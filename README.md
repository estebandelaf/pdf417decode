Pdf417decode
============

This is a fork of the project available in
<http://sourceforge.net/projects/pdf417decode>.

Usage
-----

To decode a pbm file "jac.pbm", do "./pdf417decode jac.pbm".  The file is
written to stdout.

The decoder accepts the following command line options:

 -c  output the codewords found in the image

 -d  output debug information

 -e  output the decoded information in John Lien's (jtlien@charter.net)
     pdf147_encode input format.

 -rs perform Reed-Solomon error detection and correction (note that the
     R-S error correction algorithm cannot correct certain errors
     introduced by the pbm image decoder, like the insertion of spurious
     codewords due to the same row being detected more than once).


Installation
------------

To compile the decoder, just do a "make". Assuming you have a suitable C
compiler (see Notes), the compilation will generate the executable named
pdf147decode.

In the directory "test" you'll find a few pbm test images. For each image,
there is a corresponding text file with the information used to generate the
image (in pdf417_encode input format). The file has the same name as the
image, but with "txt" extension. Most of the test images were generated with
pdf417enc 3.1, a few of them came with the pdf417enc 3.1 distribution and
one or two were found in the web.

After the program has been compiled, a "make check" can be used to test
the decoder. The command causes the program to decode all the images
found in the test directory, dumping the information to the terminal.


Notes
-----

- Your compiler needs to understand that "long long" is 64 bits (gcc does)
  in order to compile it.

- The decoder ignores Macro PDF and extended mode commands.

- The program expects the image to be in PBM (black and white) format.
  The image must be oriented horizontally, and it is processed from left to
  right and from top to bottom. So if you have an scanned image that does
  not decode into what you would expect, try then flipping it horizontally
  and/or vertically.

- The image decoder currently ignores the start and stop symbols, as well
  as the left and right row indicators. While they are not strictly necessary
  in order to decode data, their use would make the image decoder much more
  robust, and would allow the image to be decoded properly even if it is
  flipped in any direction or even rotated. They also make the error
  recovering mechanism much more accurate.

- The current Reed-Solomon implementation may fail to work properly
  for large images with error correction level greater than 4. That's
  partially related to the above limitation.

- If you use John Lien's pdf147_encode program to generate PDF417 images,
  keep in mind that the latest version 3.1 has a bug: it generates incorrect
  numeric compactions if the group of digits has more than 44 bytes.
  A simple patch can correct this, the provided pdf147_enc.patch will
  correct the offending file when you apply it on pdf417_enc.3.1 source
  directory.


History
-------

 22/12/01	Original release by Ian Goldberg. Supports only BYTE
		compactions.

 23/12/01	Modified and updated by OOO S. (ooosawaddee3@hotmail.com)
		Added a test to see if fopen() succeeded.

 07/03/04	Codeword decoding routines completely rewritten. Fixed
		a couple of bugs in dham[] table. Added support for all
		compaction modes and Reed-Solomon error correction.
		By Hector Peraza (peraza@uia.ua.ac.be).
