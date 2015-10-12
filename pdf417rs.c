#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#define NN	1024
#define PRIM	1
#define GPRIME	929
#define KK	32
#define A0	928
#define Ldec	1

#ifndef FALSE
#define FALSE	0
#endif
#ifndef TRUE
#define TRUE	1
#endif

int Alpha_to[1024];
int Index_of[1024];

static int rs_init = FALSE;


// initialize table of 3**i for syndrome calculation

void powers_init() {
    int ii;
    int power_of_3;
    int debug;

    debug = 0;
    power_of_3 = 1;
    Index_of[1] = GPRIME - 1;

    for (ii = 0; ii < GPRIME - 1; ii += 1) {
	Alpha_to[ii] = power_of_3;
	if (power_of_3 < GPRIME) {
	    if (ii != GPRIME - 1) Index_of[power_of_3] = ii;
	} else {
	    printf("Internal error: powers of 3 calculation\n");
	}

	if (debug) {
	    printf("pow = %d  ii = %d \n", power_of_3, ii);
	}
	if (debug) {
	    printf("log = %d \n", ii);
	}
	power_of_3 = (power_of_3 * 3) % GPRIME;
    }
    Index_of[0] = GPRIME - 1;
    Alpha_to[GPRIME - 1] = 1;
    Index_of[GPRIME] = A0;
}


#define modbase(x) ((x) % (GPRIME - 1))

/*
 * Performs ERRORS+ERASURES decoding of RS codes. If decoding is successful,
 * writes the codeword into data[] itself. Otherwise data[] is unaltered.
 *
 * Return number of symbols corrected, or -1 if codeword is illegal
 * or uncorrectable. If eras_pos is non-null, the detected error locations
 * are written back. NOTE! This array must be at least NN-KK elements long.
 * 
 * First "no_eras" erasures are declared by the calling program. Then, the
 * maximum # of errors correctable is t_after_eras = floor((NN-KK-no_eras)/2).
 * If the number of channel errors is not greater than "t_after_eras" the
 * transmitted codeword will be recovered. Details of algorithm can be found
 * in R. Blahut's "Theory ... of Error-Correcting Codes".
 *
 * Warning: the eras_pos[] array must not contain duplicate entries; decoder failure
 * will result. The decoder *could* check for this condition, but it would involve
 * extra time on every decoding operation.
 */

int eras_dec_rs(int data[NN], int eras_pos[NN - KK], int no_eras,
		int data_len, int synd_len)
{
    int deg_lambda, el, deg_omega;
    int i, j, r, k;
    int u, q, tmp, num1, num2, den, discr_r;
    int lambda[2048 + 1], s[2048 + 1];	/* Err+Eras Locator poly
					 * and syndrome poly */
    int b[2048 + 1], t[2048 + 1], omega[2048 + 1];
    int root[2048], reg[2048 + 1], loc[2048];
    int syn_error, count;
    int ci;
    int error_val;
    int fix_loc;
    int debug;

    if (!rs_init) powers_init();

    debug = 0;

#define DEBUG 2

#if DEBUG >= 1 && MM != 8

    /* Check for illegal input values */
    for (i = 0; i < data_len; i++)
	if (data[i] > GPRIME)
	    return -1;
#endif

    /* form the syndromes; i.e. evaluate data(x) at roots of g(x)
       namely @**(1+i)*PRIM, i = 0, ... , (NN-KK-1) */
    for (i = 1; i <= synd_len; i++) {
	s[i] = 0;//data[data_len];
    }

    for (j = 1; j <= data_len; j++) {

	if (data[data_len - j] == 0) continue;

	tmp = Index_of[data[data_len - j]];

	/*  s[i] ^= Alpha_to[modbase(tmp + (1+i-1)*j)]; */
	for (i = 1; i <= synd_len; i++) {
	    s[i] = (s[i] + Alpha_to[modbase(tmp + (i) * j)]) % GPRIME;
	}
    }

    /* Convert syndromes to index form, checking for nonzero condition */
    syn_error = 0;
    for (i = 1; i <= synd_len; i++) {
	syn_error |= s[i];
	if (debug) {
	    printf("Raw syndrome = %d i = %d \n", s[i], i);
	}
	s[i] = Index_of[s[i]];
    }

    if (!syn_error) {
	/* if syndrome is zero, data[] is a codeword and there are no
	 * errors to correct. So return data[] unmodified
	 */
	count = 0;
	printf("No errors \n");
	goto finish;
    }

    for (ci = synd_len - 1; ci >= 0; ci--) lambda[ci + 1] = 0;

    lambda[0] = 1;

    if (no_eras > 0) {
	/* Init lambda to be the erasure locator polynomial */
	lambda[1] = Alpha_to[modbase(PRIM * eras_pos[0])];
	for (i = 1; i < no_eras; i++) {
	    u = modbase(PRIM * eras_pos[i]);
	    for (j = i + 1; j > 0; j--) {
		tmp = Index_of[lambda[j - 1]];
		if (tmp != A0)
		    lambda[j] = (lambda[j] + Alpha_to[modbase(u + tmp)]) % GPRIME;
	    }
	}
#if DEBUG >= 1
	/* Test code that verifies the erasure locator polynomial just constructed
	   Needed only for decoder debugging. */

	/* find roots of the erasure location polynomial */
	for (i = 1; i <= no_eras; i++)
	    reg[i] = Index_of[lambda[i]];
	count = 0;		// usg NN = GPRIME-1?

	for (i = 1, k = data_len - Ldec; i <= data_len + synd_len;
	     i++, k = modbase(data_len + k - Ldec)) {
	    q = 1;
	    for (j = 1; j <= no_eras; j++)
		if (reg[j] != A0) {
		    reg[j] = modbase(reg[j] + j);
		    q = (q + Alpha_to[reg[j]]) % GPRIME;
		}
	    if (q != 0)
		continue;
	    /* store root and error location number indices */
	    root[count] = i;
	    loc[count] = k;
	    count++;
	}
	if (count != no_eras) {
	    // printf("\n lambda(x) is WRONG\n");
	    // count = -1;
	    //  goto finish;
	}
#if DEBUG >= 2
	printf
	    ("\n Erasure positions as determined by roots of Eras Loc Poly:\n");
	for (i = 0; i < count; i++)
	    printf("loc  = %d ", loc[i]);
	printf("\n");
#endif
#endif
    }
    for (i = 0; i < synd_len + 1; i++)
	b[i] = Index_of[lambda[i]];

    /*
     * Begin Berlekamp-Massey algorithm to determine error+erasure
     * locator polynomial
     */
    r = no_eras;
    el = no_eras;
    while (++r <= synd_len) {	/* r is the step number */
	/* Compute discrepancy at the r-th step in poly-form */
	discr_r = 0;
	for (i = 0; i < r; i++) {
	    if ((lambda[i] != 0) && (s[r - i] != A0)) {
		if (debug) {
		    printf("do add Index_of[lambda[]] = %d \n",
			   Index_of[lambda[i]]);
		}
		if (i % 2 == 1) {
		    discr_r = (discr_r + Alpha_to[modbase((Index_of[lambda[i]] + s[r - i]))]) % GPRIME;
		} else {
		    discr_r = (discr_r + GPRIME - Alpha_to[modbase((Index_of[lambda[i]] + s[r - i]))]) % GPRIME;
		}
		if (debug) {
		    printf("In loop - discr = %d i = %d r = %d lambda[i] = %d s[r-i] = %d \n",
			   discr_r, i, r, lambda[i], s[r - i]);
		}
	    }
	}
	if (debug) {
	    printf("r = %d Discrepency = %d \n", r, discr_r);
	}

	discr_r = Index_of[discr_r];	/* Index form */

	if (discr_r == A0) {
	    /* 2 lines below: B(x) <-- x*B(x) */
	    //  COPYDOWN(&b[1],b,synd_len);
	    //
	    if (debug) {
		printf("Discrepancy = A0\n");
	    }
	    for (ci = synd_len - 1; ci >= 0; ci--) b[ci + 1] = b[ci];
	    b[0] = A0;
	} else {
	    /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
	    /*  the T(x) will become the next lambda */

	    t[0] = lambda[0];
	    for (i = 0; i < synd_len; i++) {
		if (debug) {
		    printf("i = %d b[i] = %d \n", i, b[i]);
		}
		if (b[i] != A0) {

		    //  t[i+1] =  (lambda[i+1] + GPRIME -
		    //              Alpha_to[modbase(discr_r + GPRIME - 1 -  b[i])]) % GPRIME;
		    t[i + 1] = (lambda[i + 1] + Alpha_to[modbase(discr_r + b[i])]) % GPRIME;

		    if (debug) {
			printf("New t[i+1] = %d lambda[i+1] = %d b[i] = %d i = %d discr_r = %d\n",
			       t[i + 1], lambda[i + 1], b[i], i, discr_r);
		    }
		} else {
		    t[i + 1] = lambda[i + 1];
		}
		if (debug) {
		    printf("i = %d t[i+1] = %d lambda[i+1] = %d \n", i,
			   t[i + 1], lambda[i + 1]);
		}
	    }
	    el = 0;
	    if (2 * el <= r + no_eras - 1) {
		if (debug) {
		    printf("Reached the el stuff, inv  el = %d r = %d \n", el, r);
		}
		el = r + no_eras - el;
		/*
		 * 2 lines below: B(x) <-- inv(discr_r) *
		 * lambda(x)
		 */
		for (i = 0; i <= synd_len; i++) {

		    if (lambda[i] == 0) {
			b[i] = A0;
		    } else {
			b[i] = modbase(Index_of[lambda[i]] - discr_r + GPRIME - 1);
			if (debug) {
			    printf("Inverting le  b[i] = %d i = %d \n", b[i], i);
			}
		    }
		}

	    } else {
		if (debug) {
		    printf("Reached the el stuff, x mul,   el = %d r = %d \n", el, r);
		}
		/* 2 lines below: B(x) <-- x*B(x) */
		//      COPYDOWN(&b[1],b,synd_len);
		for (ci = synd_len - 1; ci >= 0; ci--) b[ci + 1] = b[ci];
		b[0] = A0;
	    }
	    //      COPY(lambda,t,synd_len+1);

	    for (ci = synd_len + 1 - 1; ci >= 0; ci--) {
		lambda[ci] = t[ci];
		if (debug) {
		    printf("ci = %d Lambda = %d \n", ci, t[ci]);
		}
	    }
	}
    }

    /* Convert lambda to index form and compute deg(lambda(x)) */
    deg_lambda = 0;
    for (i = 0; i < synd_len + 1; i++) {

	lambda[i] = Index_of[lambda[i]];

	if (lambda[i] != A0) deg_lambda = i;

	if (debug) {
	    printf("Lambda in index form = %d \n", lambda[i]);
	}

    }

    if (debug) {
	printf("Determination of deg_lambda = %d \n", deg_lambda);
    }

    /*
     * Find roots of the error+erasure locator polynomial by Chien
     * Search
     */

    for (ci = synd_len - 1; ci >= 0; ci--)
	reg[ci + 1] = lambda[ci + 1];

    count = 0;			/* Number of roots of lambda(x) */
    for (i = 1, k = data_len - 1; i <= GPRIME; i++) {
	q = 1;
	if (debug) {
	    printf(" Reg[j] = %d q = %d i = %d \n", reg[j], q, i);
	}
	for (j = deg_lambda; j > 0; j--) {

	    if (reg[j] != A0) {
		if (debug) {
		    printf("loop Reg[j] pre = %d \n", reg[j]);
		}
		reg[j] = modbase(reg[j] + j);
		//      q = modbase( q +  Alpha_to[reg[j]]);
		if (deg_lambda != 1) {
		    if (j % 2 == 0) {
			q = (q + Alpha_to[reg[j]]) % GPRIME;
		    } else {
			q = (q + GPRIME - Alpha_to[reg[j]]) % GPRIME;
		    }
		} else {
		    q = Alpha_to[reg[j]] % GPRIME;
		    if (q == 1) --q;
		}
		if (debug) {
		    printf("loop Reg[j] = %d q = %d i = %d j = %d %d = k\n",
			   reg[j], q, i, j, k);
		}
	    }
	}

	if (q == 0) {
	    /* store root (index-form) and error location number */
	    root[count] = i;

	    loc[count] = GPRIME - 1 - i;
	    if (count < synd_len) {
		count += 1;
	    } else {
		printf("Error : Error count too big = %d \n", count);
	    }

	    if (debug) {
		printf("root  = %d loc = %d \n", i, k);
	    }

	}
	if (k == 0) {
	    k = data_len - 1;
	} else {
	    k -= 1;
	}

	/* If we've already found max possible roots,
	 * abort the search to save time
	 */

	if (count == deg_lambda) break;

    }

    if (deg_lambda != count) {
	/*
	 * deg(lambda) unequal to number of roots => uncorrectable
	 * error detected
	 */

	printf("Uncorrectable error: root count = %d deg lambda = %d \n",
	       count, deg_lambda);
	count = -1;
	goto finish;
    }

    /*
     * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
     * x**(synd_len)). in index form. Also find deg(omega).
     */
    deg_omega = 0;
    for (i = 0; i < synd_len; i++) {
	tmp = 0;
	j = (deg_lambda < i) ? deg_lambda : i;
	if (debug) {
	    printf("j = %d deg_lambda = %d lambda[j] = %d \n",
		   j, deg_lambda, lambda[j]);
	}
	for (; j >= 0; j--) {
	    if ((s[i + 1 - j] != A0) && (lambda[j] != A0)) {
		if (j % 2 == 1) {
		    tmp = (tmp + GPRIME - Alpha_to[modbase(s[i + 1 - j] + lambda[j])]) % GPRIME;
		} else {

		    tmp = (tmp + Alpha_to[modbase(s[i + 1 - j] + lambda[j])]) % GPRIME;
		}
		if (debug) {
		    printf("In tmp loop  tmp = %d i = %d j = %d s[i+1-j] = %d lambda[j] = %d \n",
			   tmp, i, j, s[i + 1 - j], lambda[j]);
		}
	    }
	}

	if (tmp != 0) deg_omega = i;
	omega[i] = Index_of[tmp];
	if (debug) {
	    printf("Omega [i] = %d i = %d \n", omega[i], i);
	}

    }
    omega[synd_len] = A0;
    if (debug) {
	printf("Degree of omega = %d \n", deg_omega);
    }

    /*
     * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
     * inv(X(l))**(B0-1) and den = lambda_pr(inv(X(l))) all in poly-form
     */
    for (j = count - 1; j >= 0; j--) {
	num1 = 0;
	for (i = deg_omega; i >= 0; i--) {
	    if (omega[i] != A0) {
		//    num1  = ( num1 + Alpha_to[modbase(omega[i] + (i * root[j])]) % GPRIME;
		num1 = (num1 + Alpha_to[modbase(omega[i] + ((i + 1) * root[j]))]) % GPRIME;
		if (debug) {
		    printf("Num1 = %d i = %d omega[i] = %d root[j] = %d \n",
			   num1, i, omega[i], root[j]);
		}
	    }
	}
	//  num2 = Alpha_to[modbase(root[j] * (1 - 1) + data_len)];

	num2 = 1;
	den = 0;

	// denominator if product of all (1 - Bj Bk) for k != j
	// if count = 1, then den = 1

	den = 1;
	for (k = 0; k < count; k += 1) {
	    if (k != j) {
		tmp = (1 + GPRIME - Alpha_to[modbase(GPRIME - 1 - root[k] + root[j])]) % GPRIME;
		den = Alpha_to[modbase(Index_of[den] + Index_of[tmp])];
	    }
	}

	if (debug) {
	    printf("den = %d \n", den);
	}

	if (den == 0) {
#if DEBUG >= 1
	    printf("\n ERROR: denominator = 0\n");
#endif
	    /* Convert to dual- basis */
	    count = -1;
	    goto finish;
	}

	if (debug) {
	    printf("Index num1 = %d Index num2 = %d Index of den = %d \n",
		   Index_of[num1], Index_of[num2], Index_of[den]);
	}

	error_val = Alpha_to[modbase(Index_of[num1] + Index_of[num2] +
				     GPRIME - 1 - Index_of[den])] % GPRIME;

	if (debug) {
	    printf("error_val = %d \n", error_val);
	}

	/* Apply error to data */
	if (num1 != 0) {
	    if (loc[j] < data_len + 1) {
		fix_loc = data_len - loc[j];
		if (debug) {
		    printf("Fix loc = %d \n", fix_loc);
		}
		if (fix_loc < data_len + 1) {
		    data[fix_loc] = (data[fix_loc] + GPRIME - error_val) % GPRIME;
		}
	    }
	}
    }
  finish:
    if (debug) {
	printf("At FINISH \n");
    }
    if (eras_pos != NULL) {
	for (i = 0; i < count; i++) {
	    if (eras_pos != NULL) eras_pos[i] = loc[i];
	}
    }

    return count;
}
