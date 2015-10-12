/* Pdf417decode.c

   Pdf417decode can decode PDF417 barcode information from a pbm file.

   Original software by Ian Goldberg, 
   Modified and updated by OOO S.  (ooosawaddee3@hotmail.com)
   Version 2.0 by Hector Peraza (peraza@uia.ua.ac.be)

   Usage: To decode a pbm file "jac.pbm", do "./pdf417decode jac.pbm".  
   The file is written to stdout.

   See the accompaining Readme.txt file for information about the
   command line options.
   
   Notes:

   - Your compiler needs to understand that "long long" is 64 bits (gcc does) 
     in order to compile it.

   If you can debug or modify this software to work better than now 
   (like add error correction or modify for any operation system) 
   you can cantact OOO S.  (ooosawaddee3@hotmail.com)
   I will annouce your name in the next modify rev.

   For specification of pdf417 you can see at
      http://www.geocities.com/ooosawaddee3pdf417/ 

   History:

   22/12/01  pdf417decode          Initial Release 
   23/12/01  pdf417decode rev 1.1  Added a test to see if fopen() suceeded.
   07/03/04  pdf417decode rev 2.0  Decoding routines heavily modified.
                                   Fixed start-of-row symbol problem.
                                   Numeric and text compactions now supported.

*/


#include <stdio.h>
#include <pbm.h>


/* You may have to play with these numbers, depending on your scan quality */

#define FUZZ_THRESH (cols/40)
#define ROW_THRESH  (cols/20)

typedef unsigned int UInt32;
typedef int Int32;

/*  You'd better be using gcc, or else an Alpha, or something
 *  Your compiler needs to understand that "long long" is 64 bits (gcc does) 
 *  in order to compile it.
 */

typedef unsigned long long UInt64;

extern const unsigned dham[3][32768];

int mask[15];

int debug = 0;
int dump = 0;
int encfmt = 0;

Int32 codewords[34*90];  /* array for the extracted codewords */
Int32 erasures[34*90];   /* for the Reed-Solomon correction routine */
int numouts = 0;
int numerasures = 0;

void decode_segment(int *cw, int len, int mode);

void convert_byte(int *cw, int len, int mode);
void convert_text(int *cw, int len);
void convert_num(int *cw, int len);

extern int eras_dec_rs(int data[], int eras_pos[], int no_eras,
                       int data_len, int synd_len);
                


/*-----------------------------------------------------------------*/

/* Return -1 if wrong cluster */

static int bestham(int word, int which) {
    int i, best = which;

    if (dump) printf("%d %.4x 0x%08x (%d)\n",
                     which, word, dham[which][word], dham[which][word] & 0xffff);

    for (i = 0; i < 3; ++i) {
	if ((dham[i][word] & 0xff000000) < (dham[best][word] & 0xff000000)) {
	    best = i;
	}
    }
    if (best != which) return -1;

    return dham[which][word];
}


static void add_codeword(int word) {
    static int len = 0;
    static int sorow = 0;  /* start of row */
    static int skip = 0;

    if (skip) {
	--skip;
	return;
    }

    if (numouts == 0) {
	len = word;  /* not really used, but stored as codewords[0] */
    }

    if ((word & 0xffffff) == 0x030000) {  /* start sequence */
    //if ((word & 0xff00000) == 0x2000000) {  /* start sequence */
	skip = 1;
	sorow = numouts;
	return;
    }

    if ((word & 0xffffff) == 0x030001) {  /* stop sequence? */
	--numouts;
	return;
    }

    if (word == -3) {
	/* Rewind to beginning of row */
	numouts = sorow;
	return;
    }

    if ((word & 0xffff) == 0xffff) {
	word = 0;
	erasures[numerasures++] = numouts;
    }

    /* store codeword */
    codewords[numouts++] = word & 0xffff;
}


/*
 *  This routine does all the decoding. Individual compaction segments
 *  are extracted from the codeword array, then for each segment the
 *  decode_segment() function is called, which in turn invokes one of
 *  the convert_*() functions which finally decode the data into a
 *  human readable format.
 */

static void decode_codewords() {
    int   i, cw, len, slen, mode, shift;
    Int32 segment[34*90];  /* single compaction segment to be decoded */

    if (numouts == 0) return;
    
    len = codewords[0];
    if (len == 0) return;

    slen = 0;

    mode = 900;    /* default mode is Text Compaction */
    shift = mode;
    
    for (i = 1; i < len; ++i) {
        cw = codewords[i];
        if (cw >= 900) {
        
            if (slen > 0) decode_segment(segment, slen, mode);
            slen = 0;
            
            switch (cw) {
            case 900:  /* mode latch to Text Compaction mode */
            case 901:  /* mode latch to Byte Compaction */
            case 902:  /* mode latch to Numeric Compaction */
		mode = shift = cw;
		break;
		
            case 913:  /* mode shift to Byte Compaction */
            	shift = cw;
            	break;
            
            case 921:  /* reader initialization */
            	break;
            	
            case 922:  /* terminator codeword for Macro PDF control block */
                break;
                
            case 923:  /* Begin of optional fields in the Macro PDF control block */
                break;
                
            case 924:  /* mode latch to Byte Compaction (num of encoded bytes
                          is an integer multiple of 6) */
                mode = shift = cw;
                break;
                
            case 925:  /* identifier for a user defined Extended Channel
                          Interpretation (ECI) */
            case 926:  /* identifier for a general purpose ECI format */
            case 927:  /* identifier for an ECI of a character set or code page */
                break;
                
            case 928:  /* Begin of a Macro PDF Control Block */
                break;
                
            default:
                fprintf(stderr, "Unknown mode %d\n", cw);
                break;
            }
            continue;
	}

	segment[slen++] = cw;
	
	if (shift != mode) {
	    if (slen > 0) decode_segment(segment, slen, shift);
	    slen = 0;
	    shift = mode;
	}
    }

    if (slen > 0) decode_segment(segment, slen, mode);
}


void decode_segment(int *cw, int len, int mode) {

    switch (mode) {
    case 900:
	convert_text(cw, len);
	break;
	        
    case 901:
    case 913:
    case 924:
        convert_byte(cw, len, mode);
    	break;
	    	
    case 902:
        convert_num(cw, len);
    	break;
    }
}


void convert_byte(int *cw, int len, int mode) {
    UInt64 codeval;
    int i, j;
    unsigned char b[6];

    if (debug > 1) printf("convert_byte: %d codewords (mode = %d)\n", len, mode);

    if (encfmt) printf("BC \"");

    /* 6 bytes are encoded in a group of 5 codewords */
    for ( ; (mode == 901) ? (len > 5) : (len >= 5); len -= 5) {
        codeval = 0;

        /* convert from base 900 to base 256 */

	for (i = 0; i < 5; ++i) {
	    codeval *= 900;
	    codeval += *cw++;
	}
    
	if (debug > 1) printf("codeval = %Lx, giving ", codeval);

	for (j = 0; j < 6; ++j) {
	    b[5-j] = codeval % 256;
	    if (debug > 1) printf("[%02x] ", b[5-j]);
	    codeval >>= 8;
	}
	if (debug > 1) printf("\n");
	
	if (encfmt) {
	    for (j = 0; j < 6; ++j) printf("%02X", b[j]);
	} else {
	    for (j = 0; j < 6; ++j) printf("%c", b[j]);
	}
    }
    
    /* remaining codewords, if any, are encoded 1 byte per codeword */
    if (len > 0) {
	if (debug > 1) printf("remaining %d codewords: ", len);
	for (j = 0; j < len; ++j) {
	    b[j] = *cw++;
	    if (debug > 1) printf("[%02x] ", b[j]);
        }
	if (debug > 1) printf("\n");
	
	if (encfmt) {
	    for (j = 0; j < len; ++j) printf("%02X", b[j]);
	} else {
	    for (j = 0; j < len; ++j) printf("%c", b[j]);
	}
	
	codeval = 0;
	i = 0;
    }

    if (encfmt) printf("\"\n");

    fflush(stdout);
}


void convert_text(int *cw, int len) {
    int mode, shift, enc;
    int i, j, c[2], cout;

    static char txt_upper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ    ";
    static char txt_lower[] = "abcdefghijklmnopqrstuvwxyz    ";
    static char txt_mixed[] = "0123456789&\r\t,:#-.$/+%*=^     ";
    static char txt_punct[] = ";<>@[\\]_`~!\r\t,:\n-.$/\"|*()?{}' ";

    if (debug > 1) printf("convert_text: %d codewords\n", len);

    mode = shift = 0;

    if (encfmt) printf("TC \"");
    
    for (i = 0; i < len; ++i) {

	c[0] = *cw / 30;
	c[1] = *cw % 30;
	++cw;

	for (j = 0; j < 2; ++j) {

	    if (debug > 1) printf(" (%d)-", c[j]);

	    enc = mode;
            if (mode != shift) { enc = shift; shift = mode; }

	    switch (enc) {
	    case 0:  /* uppercase */
	        if (c[j] == 27) { mode = shift = 1; continue; }  /* lower latch */
	        if (c[j] == 28) { mode = shift = 2; continue; }  /* mixed latch */
	        if (c[j] == 29) { shift = 3; continue; }         /* punct shift */
	        cout = txt_upper[c[j]];
	        break;

	    case 1:  /* lowercase */
	        if (c[j] == 27) { shift = 0; continue; }         /* upper shift */
	        if (c[j] == 28) { mode = shift = 2; continue; }  /* mixed latch */
	        if (c[j] == 29) { shift = 3; continue; }         /* punct shift */
	        cout = txt_lower[c[j]];
	        break;

	    case 2:  /* mixed (numeric) */
	        if (c[j] == 25) { mode = shift = 3; continue; }  /* punct latch */
	        if (c[j] == 27) { mode = shift = 1; continue; }  /* lower latch */
	        if (c[j] == 28) { mode = shift = 0; continue; }  /* upper latch */
	        if (c[j] == 29) { shift = 3; continue; }         /* punct shift */
	        cout = txt_mixed[c[j]];
	        break;

	    case 3:  /* punctuation */
	        if (c[j] == 29) { mode = shift = 0; continue; }  /* upper latch */
	        cout = txt_punct[c[j]];
	        break;
	    }
    
	    printf("%c", cout);
        }

    }

    if (encfmt) printf("\"\n");
    fflush(stdout);
}


#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

/*
 *  BCD arithmetic is used in this routine in order to simplify
 *  large-precision number manipulations.
 */

void convert_num(int *cw, int len) {
    int n_bcd[45], cw_bcd[3];
    int i, j, n, res, carry, start;
    
    if (debug > 1) printf("convert_num: %d codewords\n", len);

    for ( ; len > 0; len -= 15) {

        /* clear the accumulator */
        for (i = 0; i < 45; ++i) n_bcd[i] = 0;
        
        for (i = 0; i < MIN(len, 15); ++i) {

	    /* convert codeword to BCD */
	    n = *cw++;
	    cw_bcd[0] = n % 10;
	    n /= 10;
	    cw_bcd[1] = n % 10;
	    n /= 10;
	    cw_bcd[2] = n;

	    /* multiply accumulator by 900 (100 * 9) */
            if (i > 0) {
                carry = 0;
        
                /* multiply by 9 */
                for (j = 0; j < 45; ++j) {
                    res = n_bcd[j] * 9 + carry;
                    n_bcd[j] = res % 10;
                    carry = res / 10;
                }
        
                /* multiply by 100 */
                for (j = 44; j >= 2; --j) {
                    n_bcd[j] = n_bcd[j-2];
                }
                n_bcd[1] = 0;
                n_bcd[0] = 0;
            }

            /* then add the BCD codeword */
            carry = 0;
    
            for (j = 0; j < 3; ++j) {
                res = n_bcd[j] + cw_bcd[j] + carry;
                n_bcd[j] = res % 10;
                carry = res / 10;
            }
            for ( ; j < 45; ++j) {
                res = n_bcd[j] + carry;
                n_bcd[j] = res % 10;
                carry = res / 10;
            }

        }

        start = 0;
        
        if (encfmt) printf("NC \"");
        for (j = 0; j < 45; ++j) {
            if (start) {
                printf("%c", n_bcd[44-j] + '0');
            } else if (n_bcd[44-j] == 1) {
                start = 1;
            } else if (n_bcd[44-j] != 0) {
                printf("<invalid>");
                break;
            }
        }
        if (encfmt) printf("\"\n");
    }

}


/* this routine extracts the codewords from a single pixel row from the image */

int processrow(int cols, int rownum, int num, double *cumbits) {
    int firstblack = 0;
    int scale;
    int j;
    int nchange;
    int cumchange[cols];
    double thresh = 0.5 * num;

    while ((firstblack < cols) && (cumbits[firstblack] < thresh))
	++firstblack;

    if (firstblack+1 >= cols) return 0;
    nchange = 0;
    cumchange[nchange++] = 0;

    for (j = firstblack+1; j < cols; ++j) {
	if ((cumbits[j] < thresh) != (cumbits[j-1] < thresh)) {
	    if (nchange > 1 && 
		(j - firstblack - cumchange[nchange-1]) * 15 < cumchange[1]) {
		/* Spurious change */
		--nchange;
	    } else {
		cumchange[nchange++] = j - firstblack;
	    }
	}
    }

    if (debug > 1) {
        for (j = 0; j < nchange; ++j) {
	    printf("%3d ", cumchange[j]);
        }
        printf("\n");
    }
    
    if (nchange < 8) return 0;

    for (j = 0; j < nchange-8; j += 8) {
	int k;
	int word = 0;
	unsigned cw;

	scale = cumchange[j+8] - cumchange[j];
	for (k = 0; k < 8; k += 2) {
	    int s, e, l;

	    s = (int) (17.0 * (double) (cumchange[j+k] - cumchange[j]) / (double) scale + 0.5);
	    e = (int) (17.0 * (double) (cumchange[j+k+1] - cumchange[j]) / (double) scale + 0.5);
	    /* We know we always start with 1, end with 0 */
	    if (s < 1) s = 1;
	    if (e > 16) e = 16;
	    for (l = s; l < e; ++l) word |= mask[l-1];
	}
	word >>= 1;
	cw = bestham(word, rownum % 3);
	add_codeword(cw);
    }

    return 1;
}


int main(int argc, char **argv) {
    FILE *pf;
    int rows, cols;
    bit **bits;
    double *cumbits;
    int i, j;
    int ready, num, ecc = 0;
    int rownum;
    char *myname = argv[0];

    for (i = 0; i < 15; ++i) mask[i] = 1 << (15-i);
    
    for ( ; argc > 1; --argc, ++argv) {
        if (strcmp(argv[1], "-d") == 0)
            ++debug;
        else if (strcmp(argv[1], "-c") == 0)
            dump = 1;
        else if (strcmp(argv[1], "-e") == 0)
            encfmt = 1;
        else if (strcmp(argv[1], "-rs") == 0)
            ecc = 1;
        else
            break;
    }
    pbm_init(&argc, argv);
    
    if (argc > 1) {
      pf = fopen(argv[1], "r");
    } else {
      fprintf(stderr, "usage: %s [-d] [-c] [-e] [-rs] file\n", myname);
      exit(1);
    }

    if (pf == NULL) {
      fprintf(stderr, "%s: could not open file: %s\n", myname, argv[1]);
      exit(1);
    }
 
    bits = pbm_readpbm(pf, &cols, &rows);
    fclose(pf);

    cumbits = malloc(cols * sizeof(double));

    ready = 1;
    num = 0;
    rownum = 0;

    for (j = 0; j < cols; ++j) cumbits[j] = 0.0;

    for (i = 1; i < rows; ++i) {
	int d = 0;
	for (j = 0; j < cols; ++j) {
	    if (bits[i][j] != bits[i-1][j]) ++d;
	}
	if (d < FUZZ_THRESH) {
	    if (ready == 1) {
		num = 0;
		for (j = 0; j < cols; ++j) cumbits[j] = 0.0;
		ready = 2;
	    }
	    for (j = 0; j < cols; ++j) {
		if (bits[i][j] == PBM_BLACK) cumbits[j] += 1.0;
	    }
	    ++num;
	} else if (d > ROW_THRESH) {
	    if (ready == 2) {
		if (processrow(cols, rownum, num, cumbits)) ++rownum;
		ready = 1;
	    }
	}
    }
    if (ready == 2) if (processrow(cols, rownum, num, cumbits)) ++rownum;

    if (ecc) {
	printf("Total codewords = %d (%d data, %d ECC)\n",
	       numouts, codewords[0], numouts - codewords[0]);

	//num = eras_dec_rs(codewords, erasures, numerasures, numouts, numouts - codewords[0]);
	num = eras_dec_rs(codewords, NULL, 0, numouts, numouts - codewords[0]);
	if (num < 0)
	    printf("Errors detected, but data could not be corrected\n");
	else if (num > 0)
	    printf("%d codewords corrected\n\n", num);
    }
    
    decode_codewords();

    pbm_freearray(bits, rows);
    return 0;
}
