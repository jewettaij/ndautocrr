[![GitHub repo size](https://img.shields.io/github/repo-size/jewettaij/ndautocrr)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()



ndautocrr
===========

This program calculates the autocorrelation of
a series of numbers or vectors (**x**(i)) from a text file
(which may contain one or more data sets).
It prints the autorcorrelation function

*C(j)* = ⟨(**x**(i)-⟨**x**⟩)⋅(**x**(i+j)-⟨**x**⟩)⟩

as a function of j, to the standard output, where ⟨⟩ denotes the
[average](https://en.wikipedia.org/wiki/Average#Arithmetic_mean),
and ⋅ denotes multiplication *(either scalar multiplication or
the dot product)*.




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
(If the file has more than one column (as in the example below)
the **x** are assumed to be vectors.)
to multiply different data entries together.
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
          < inputlist.txt > corrfunc.txt
```


### Notes


1. *L* is the size of the domain for the autocorrelation function.
The autocorrelation function is truncated beyond this point.
You can manually specify this by using the "-L" argument.  (See below.)
If left unspecified, **by default, the autocorrelation function is truncated
when it decays below a threshold value (which is 1/e by default)
relative to the peak, *C(0)*.**
(See the discussion of the "-t" argument below.)

2. The data entries **x**(i) can be scalar values or vectors.
When they are vectors, the
[dot product](https://en.wikipedia.org/wiki/Dot_product)
is used. (See [inner_product.h](./src/inner_product.h).)

3. For convenience, a crude estimate of the correlation length
(or correlation time) is printed to the standard error.
  * By default, the correlation length is estimated by
finding the first location *j* where *C(j)/C(0)* crosses *1/e*.
  * If the threshold is adjusted (using the *-t* argument), then the correlation
length is estimated by fitting the curve to an exponentially decaying function,
considering only the point where it crosses the user-specified threshold.
(So if the threshold is set to *1/e^2*, and *C(j)/C(0)=1/e^2*,
 then the correlation length will be reported as *j/2*.)
  * If the -L argument (*L*) is specified, the correlation length is estimated
by summing *Σ_j C(j)* from *j=0* to *j=L*.
*(Graphing the autocorrelation function is always a good way to choose an
appropriate -L parameter, especially if oscillations in the data are present.
Alternatively, you can determine L by running the program once with the "-t"
argument and a suitable threshold, counting the number of lines of output,
and then running the program again, setting the -L argument to that number.)*
  * If both -L and -t are specified, the threshold method is preferentially used
*(unless C(j)/C(0) fails to decay enough before j reaches L).*

There are much more robust methods for estimating the correlation length,
but these methods are the simplest.


## Optional arguments

The autocorrelation function *C(j)* typically decays as j increases.
In addition to this decay of signal strength, there are often
large fluctuations in the autocorrelation function *C(j)* for large
values of *j*, due to increasingly poor sampling at large *j*.
At some point, one must make a decision where to truncate *C(j)*.
You can either do this by manually specifying the maximum value
of *j* (using the "domainwidth" argument to directly specify *L*),
and/or by using a threshold cutoff ("-t").


### -L domainwidth

You can manually specify *L*, the width of the domain of *C(j)*,
by supplying an integer as one of the command line arguments.
This integer determines where the correlation function *C(j)* is truncated.
(Values of *C(j)* for *j>L* will not be computed.)
By default, *L* is determined by the point when *C(j)*
decays to *1/e* of its original value.
But you can override this choice and force *L* 
to be any number in the range from *0* to *N-1*.
Note that if "-t" argument is unspecified, then the correlation length
will be estimated by computing the sum, *Σ_j C(j)* from j=0 to j=L,
whenever "-L" is used.  *(Due to poor sampling at large L, this method for
estimating the correlation length can be very sensitive to the choice of L,
so it is not used by default.)*


### -t threshold

Again *C(j)* will typically decay as *j* increases.
By default, ndautocrr to halt when *C(j)/C(0)* decays below some
halt once *C(j)/C(0)* < *threshold*, where *threshold* is specified by the user.
If unspecified, by default this threshold is *1/e*,
The choice of threshold may also effect the estimate of the correlation length.
*(Unless the "-L" argument is used, the correlation length is estimated
by fitting of the autocorrelation function to a decaying exponential
at the location where it crosses the threshold.)*


### -ave

This means subtract the average (⟨**x**⟩) before calculating the
autocorrelation function *C(j)*.  (This is the default behavior.)


### -avezero

The "-avezero" argument allows you to force ⟨**x**⟩=0 when computing *C(j)*.
This can be useful if you are calculating the autocorrelation function
of a sequence of data **x** whose average value you know should be zero.
In most cases, the average **x** in your files will not to be
zero due to limited sampling or bias in the input file.
The "-avezero" argument can be useful in these situations.


### -p

If you pass the "-p" flag, then periodic boundary conditions
are applied to the data.
That means *(in the simple case, for a single data set of length N)*
the autocorrelation function, *C(j)*, is computed this way:

<img src="http://latex.codecogs.com/gif.latex?\large&space;C(j)=\frac{1}{N}\sum_{i=1}^{N}(\mathbf{x}(i)-\langle \mathbf{x}\rangle)\cdot(\mathbf{x}((i+j)\%N)-\langle \mathbf{x}\rangle)"/>

Here, "(i+j) % N" denotes the remainder after division by N.
(N is the number of entries in that data set.)
*(As before, if multiple data sets are used, 
then we sum over all of them when performing the average.)*


### -nsum
Report an additional column in the output file (after *C(j)*)
storing the number of entries in the sum that was used to calculate *C(j)*.
*(If there is only one data set, this is N-j, where "N" is the size of that
data set.  If multiple data sets were used, then this equals
Σ_k max(0,N_k-j)), where Σ_k denotes the sum over k, and
N_k denotes the number of entries in the kth data set.)*


### -rms
Report an additional column in the output file (after *C(j)*) storing the
root-mean-squared value of (**x**(i)-⟨**x**⟩)⋅(**x**(i+j)-⟨**x**⟩),
(considering various values of *i*, for a given *j*).
*(It is not clear to me whether this quantity is ever useful.)*


## Compilation

## Linux and Apple macOS:

```
cd src
source setup_gcc.sh   #(Apple users may prefer using "setup_clang.sh" instead.)
make
```

*(Note:  If you are not using the bash shell,
enter "bash" into the terminal beforehand.)*

*(Note: Apple users can install Xcode, which includes the clang compiler by default.  Alternatively, brew can be used to install a wide range of compilers and build tools.)*

## Windows 10:

Install HyperV (with linux), or the Windows Subsystem for Linux (WSL) and run

```
sudo apt-get install build-essential
```

and then follow the instructions above.
(Older windows users can install Cygwin or MinGW, or linux via virtualbox.)

## License

ndautocrr is available under the terms of the [MIT license](LICENSE.md).
