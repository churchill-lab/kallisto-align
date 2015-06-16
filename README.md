kallisto-export
===================
kallisto-export is a utility that makes use of [kallisto](http://pachterlab.github.io/kallisto/), but outputs the pseudo-alignments in a different [binary format](#binaryformat) for use with other tools such as [emase](https://github.com/churchill-lab/emase).

Requirements
-------------

: CMake version >= 2.8.12
: HDF5 C library
: zlib


Building
-------------

To build **kallisto-export** you simply clone and compile.

```
> git clone https://github.com/churchill-lab/kallisto-export.git
> cd kallisto-export
> cmake .
> make
```

>**Note:** [kallisto](http://pachterlab.github.io/kallisto/) will be downloaded and compiled for you and be located at: external/src/kallisto-build/src/kallisto

You can now either keep the [kallisto](http://pachterlab.github.io/kallisto/) and **kallisto-export** binaries here or move them to a more suitable location.


Using
-------

To run **kallisto-export** you need a [kallisto](http://pachterlab.github.io/kallisto/) index file and a FASTQ file.

A simplified example of building a [kallisto](http://pachterlab.github.io/kallisto/) index is:

```
> kallisto index -i example.idx exampla.fa
```

> Please refer to the [kallisto](http://pachterlab.github.io/kallisto/) documentation for all options.

Once the index file is created, you can run **kallisto-export**.

```
kallisto-export [OPTION...] - pseudo-align to kallisto-export format

      --help              Print help
  -l, --load              view the binary file
      --reads             create Reads version binary format, else
                          Equivalence Class
  -f, --file FASTQ FILE   Input Fastq File
  -i, --index INDEX FILE  Input Index File
  -b, --bin arg           kallisto-export output file
```

Output
--------

**kallisto-export** pseudo-aligns the reads and outputs the format into the following binary format.

name             | type|note
-----------------|-------------------|---------------------
format          |integer| 0 for reads, 1 for equivalence class
\#targets        |integer |
\#characters in target name|integer|repeated \#targets times
target name|string|repeated \#targets times
\#haplotypes|integer|
\#characters in haplotype name|integer|repeated \#haplotypes times
haplotype name|string|repeated \#haplotypes times

>Depending upon the format value, one of two following formats is appended.

**format = 0**

name             | type|note
-----------------|-------------------|---------------------
\#reads        |integer |
\#characters in read name|integer|repeated \#reads times
read name|string|repeated \#reads times
\#pseudo-alignments|integer|
read index|integer|repeated \#pseudo-alignments times
target index|integer|repeated \#pseudo-alignments times
bitwise flag|integer|repeated \#pseudo-alignments times

**format = 1**

name             | type|note
-----------------|-------------------|---------------------
\#equivalence classes        |integer |
equivalence class count|integer|1 entry per \#equivalence classes
\#pseudo-alignments|integer|
equivalence class index|integer|repeated \#pseudo-alignments times
target index|integer|repeated \#pseudo-alignments times
bitwise flag|integer|repeated \#pseudo-alignments times

>bitwise flag is an integer representing the number of haplotypes that are "turned on" for a target.

