[![GitHub repo size](https://img.shields.io/github/repo-size/jewettaij/ndautocrr)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()



ndautocrr
===========

This program calculates the auto-correlation of
a series of numbers (*x(i)*) from a text file
(which may contain one or more lists of scalars or vectors).
It prints the correlation function,
C(j) = ⟨(x(i)-⟨x⟩)⋅(x(i+j)-⟨x⟩)⟩,
as a function of j, to the standard output, where ⟨⟩ denotes the
[average](https://en.wikipedia.org/wiki/Average#Arithmetic_mean).



### Output file format:
By default, this information is printed to the standard output (terminal)
in two-column ascii format.

```
0 C(0)
1 C(1)
2 C(2)
:   :
L C(L)
```
The domain-length, *L*, is chosen automatically, but can be specified
using the "-L" and "-t" arguments.

*This was originally a crude program written to analyze polymer simulation
trajectory files.  I find myself referring to it frequently in some of the
[moltemplate](https://github.com/jewettaij/moltemplate)
examples, so I was forced to make this ugly code available.*
*This program could probably be implemented
using only a couple lines in python.*

Multiprocessor support was implemented using
[OpenMP.](https://en.wikipedia.org/wiki/OpenMP)


### Input file format:

The text files read by this program may contain multiple sets of data
(delimited by white space, see example below).
The input file should be a text file containing 1 data entry per line,
and blank lines separating different data sets.
The example below has 2 data sets with 3 columns per entry
(but any number of columns is allowed).
```
x1 y1 z1     #<--1st dataset
x2 y2 z2
:  :  :
xN1 yN1 zN1

x1 y1 z1     #<--2nd dataset
x2 y2 z2
:  :  :
xN2 yN2 zN2
```
*When the input file contains multiple data sets, the sum used when computing averages weights each entry (each line in the file) equally.*
*Comments in the input stream (following the \# character) are ignored.*

## Usage:

```
ndautocrr [-L domainwidth] [-p] [-t threshold] [-ave,-avezero] \
          < inputlist.dat > corrfunc.dat
```



### Notes


1. *L* is the size of the domain for the auto-correlation function.
(By default, L is determined automatically.  See the discussion of the
"-t" argument below.)  However you can manually override this choice
 by using the "-L" argument.  See below.)

2. The data entries *x(i)* can be scalar values or vectors.
When they are vectors, the
[dot-product](https://en.wikipedia.org/wiki/Dot_product)
is used.

3. For convenience, the correlation length (also called the correlation time, *(sum_j C(j)/C(0))* is printed to the standard error:  


## Optional arguments

The correlation function C(j) typically decays as j increases.
In addition to this decay of signal strength, there are often
large fluctuations in the correlation function C(j) for large
values of j, due to increasingly poor sampling at large j.
At some point, one must make a decision where to truncate C(j).
You can either do this by manually specifying the maximum value
of j (using the "domainwidth" argument to directly specify L),
and/or by using a threshold cutoff ("-t").


### -L domainwidth

You can manually specify L, the width of the domain of C(j),
by supplying an integer as one of the command line arguments.
This integer determines the number of lines in the output file
(See "L" in output file format above.)  L also determines the
number of terms that will be used to calculate the correlation
length/time according to the formula shown above.
By default, L is ⌊N/2⌋, where N is the number of entries in the data set
(where ⌊⌋ denotes the 
[floor function](https://en.wikipedia.org/wiki/Floor_and_ceiling_functions)).
But you can override this choice and force L 
to be any number in the range from 0 to N-1.
A low value of L can speed speed up the program 
since the running time is O(N\*L) (for a single data set).
(Note:  If L approaches or exceeds N, a warning message will be generated.)


### -t threshold

Again C(j) will typically decay as j increases.
Alternatively, you can tell ndautocrr to halt when
C(j) decays below some threshold that the user has specified.
the width of the correlation function by using a threshold cutoff.
(By default this threshold is 0, since in many cases, negative
correlations, ie C(j) < 0, are not physically meaningful and
signify that noise has exceeded the magnitude of the signal.)
The "-t" argument should be followed by a floating point number
("threshold"), a number between -1.0 and 1.0. This number is
to be expressed as a _fraction_.  Calculation of the
autocorrelation-function will halt once C(j)/C(0) < threshold.
Of course, the subsequent calculation of the correlation
length/time only considers terms C(j) which exceed this threshold.

**Note: Use "-t -1" if you want to disable thresholding behavior.**
In that case, the default value of L is ⌊N/2⌋, where N is the
number of entries in the data set (where ⌊⌋ denotes the 
[floor function](https://en.wikipedia.org/wiki/Floor_and_ceiling_functions)).


### -ave

This means subtract the average before calculating the
correlation function C(j).  (This is the default behavior.)


### -avezero

If this argument is passed, then DO NOT subtract the average before calculating
the covariance.  This can be useful if you are calculating the correlation
length of a quantity whose average value you know to be zero (such as the
cosign of the angle between two vectors of random orientation), but whose
average turns out not to be zero due to statistical fluctuations.  In this
case you can force the average to be zero using the "-avezero"
argument.


### -p

If you pass the "-p" flag, then periodic boundary conditions
are applied to the data.
That means *(in the simple case, for a single data set of length N)*
the correlation function, C(j), is computed this way:
```
                                           __N__
                                       1   \
 C(j) = ⟨ (x(i)-⟨x⟩)⋅(x(i+j)-⟨x⟩) ⟩ = ---   >    (x(i)-⟨x⟩)⋅(x((i+j % N)-⟨x⟩)
                                       N   /____
                                            i=1
```
Here, "(i+j) % N" denotes the remainder after division by N.
(N is the number of entries in that data set.)
*(As before, if multiple data sets are used, 
  then we sum over all of them when performing the average.)*


### -nsum
Report an additional column in the output file (after C(j)) storing the number of entries in the sum that was used to calculate C(j).  *(If there is only one data set, this is N-j, where "N" is the size of that data set.  If multiple data sets were used, then this is is sum_k max(0,N_k-j)), where N_k denotes the number of entries in the kth data set.)*


### -rms
Report an additional column in the output file (after C(j)) storing the root-mean-squared value of (x(i)-⟨x⟩)⋅(x(i+j)-⟨x⟩), (considering various values of i, for a given j).*(It is not clear to me whether this quantity is ever useful.)*


## Compilation

## Linux and Apple macOS:

```
cd src
source setup_gcc.sh
make
```

*(Note:  If you are not using the bash shell,
enter "bash" into the terminal beforehand.)*

*(Note: Apple users will need to install the gcc compiler
and other build tools using Xcode or brew.)*

*(Note: If you receive an error regarding "omp" or "OpenMP", then use
"setup_gcc_serial.sh" instead.  This may be necessary for apple users.)*

## Windows 10:

Install HyperV (with linux), or the Windows Subsystem for Linux (WSL) and run

```
sudo apt-get install build-essential
```

and then follow the instructions above.
(Older windows users can install Cygwin or MinGW, or linux via virtualbox.)

## License

ndautocrr is available under the terms of the [MIT license](LICENSE.md).
