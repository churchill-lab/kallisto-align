**kallisto-align** runs [kallisto](http://pachterlab.github.io/kallisto/), processes its internal data object into an alignment incidence matrix, and exports the matrix in a pre-defined [binary format](#output) for use with other tools such as [emase](https://github.com/churchill-lab/emase) or [emase-zero](https://churchill-lab.github.io/emase-zero). It is a companion utility for our bigger project that aims to provide a suite of tools to compare, contrast, or combine the results of different alignment strategies as envisioned in our [alntools](https://churchill-lab.github.io/alntools).

* Free software: GPLv3 license

Requirements
------------

[CMake version >= 2.8.12](https://cmake.org/download/)<br />
[HDF5 C library](https://support.hdfgroup.org/HDF5/release/obtainsrc.html)<br />
[zlib](http://zlib.net)


Installation
------------

To build ```kallisto-align``` you simply clone and compile.

```
$ git clone https://github.com/churchill-lab/kallisto-align.git
$ cd kallisto-align
$ cmake .
$ make
```

>**Note:** [kallisto](http://pachterlab.github.io/kallisto/) will be downloaded and compiled for you and be located at: ```external/src/kallisto-build/src/kallisto```.

You can now either keep the ```kallisto``` and ```kallisto-align``` binaries where they are or move them to a more suitable location, e.g. ```/usr/local/bin```.


Usage
-----

We first need to build a [kallisto](http://pachterlab.github.io/kallisto/) index file. For example, given our target sequences are in a fasta file ```transcriptome.fa```:

```
$ kallisto index -i transcriptome.idx transcriptome.fa
```

> Please refer to the [kallisto manual](https://pachterlab.github.io/kallisto/manual) for more options.

Once the index file is created, you can run ```kallisto-align```:

```
kallisto-align [OPTIONS...]

Options:

  -f, --file READ FILE    Input read file name
  -i, --index INDEX FILE  Input index file name
  -b, --bin OUTPUT FILE   Output file name
  --reads                 Toggle read-level binary format (Format=0), else Equivalence Classes
  -l, --load              View the binary file
  --help                  Print this message and exit
```

For example,

```
$ kallisto-align -f my_sample.fastq -i transcriptome.idx -b my_sample.bin
```

The binary file can be converted into [emase](https://github.com/churchill-lab/emase) format with [alntools](https://churchill-lab.github.io/alntools/#ec2emase). It can also be used directly with [emase-zero](https://churchill-lab.github.io/emase-zero).


Output
------

![](assets/images/emase_binary_format.jpg "EMASE binary format")
