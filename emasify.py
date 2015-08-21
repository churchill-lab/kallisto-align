#!/usr/bin/env python
import sys
import time
import argparse
import numpy as np
from collections import OrderedDict
from emase import AlignmentPropertyMatrix as APM
from progressbar import ETA, Percentage, ProgressBar, Bar


def simple_from_one(c, size):
    ret = [0]*size
    for i in xrange(0, size):
        if (c & (1 << i)) != 0:
            ret[i] = 1
    return ret


def emasify(binary_file_name, emase_file_name):
    """

    :param binary_file_name:
    :param emase_file_name:
    :return:
    """

    if not binary_file_name:
        raise ValueError("empty file name, cannot load")

    print "Binary File: {0}".format(binary_file_name)

    f = open(binary_file_name, 'rb')

    file_version = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]

    if file_version == 0:
        print "Version: 0, Reads"
    elif file_version == 1:
        print "Version: 1, Equivalence Class"
    else:
        print "Unknown version, exiting"

    # TARGETS

    target_ids = []
    targets = OrderedDict()

    num_targets = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    print "Target Count: {0}".format(num_targets)

    for i in xrange(0, num_targets):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        target = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        targets[target] = i
        target_ids.append(target)

    # HAPLOTYPES

    haplotype_ids = []
    haplotypes = OrderedDict()

    num_haplotypes = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    print "Haplotype Count: {0}".format(num_haplotypes)

    for i in xrange(0, num_haplotypes):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        haplotype = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        haplotypes[haplotype] = i
        haplotype_ids.append(haplotype)

    if file_version == 0:

        # READS

        read_ids = []
        reads = OrderedDict()

        num_reads = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        print "Read Count: {0}".format(num_reads)

        for i in xrange(0, num_reads):
            str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
            read_id = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
            reads[read_id] = i
            read_ids.append(read_id)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        print "Alignment Count: {0}".format(num_alignments)

        alignments = np.fromfile(f, dtype = np.dtype('i'), count=num_alignments*3)

        print 'Creating APM...'
        new_shape = (len(target_ids), len(haplotypes), len(read_ids))
        aln_mat_kallisto = APM(shape=new_shape, haplotype_names=haplotype_ids, locus_names=target_ids, read_names=read_ids)

        print 'Parsing alignments...'
        widgets = [Bar('>'), ' ', ETA(), ' ', Percentage()]
        pbar = ProgressBar(widgets=widgets, maxval=num_alignments*3).start()
        counter = 0

        for i in xrange(0, num_alignments*3, 3):
            rid = alignments[i]
            lid = alignments[i+1]
            temp_bits = alignments[i+2]

            counter += 1
            pbar.update(i)
            if temp_bits == 0:
                continue

            bits = simple_from_one(temp_bits, num_haplotypes)
            for hid, b in enumerate(bits):
                if b:
                    aln_mat_kallisto.set_value(lid, hid, rid, 1)

        pbar.finish()
        print "Finalizing..."
        aln_mat_kallisto.finalize()
        aln_mat_kallisto.save(emase_file_name, title='KALLISTOALIGN')

        print "DONE"
    else:

        # EQUIVALENCE CLASSES

        num_ec = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        print "Equivalance Class Count: {0}".format(num_ec)

        ec_ids = [x for x in xrange(0, num_ec)]
        counts = np.fromfile(f, dtype=np.dtype('i'), count=num_ec)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        print "Alignment Count: {0}".format(num_alignments)

        alignments = np.fromfile(f, dtype = np.dtype('i'), count=num_alignments*3)

        print 'Creating APM...'
        new_shape = (len(target_ids), len(haplotypes), len(counts))
        aln_mat_kallisto = APM(shape=new_shape, haplotype_names=haplotype_ids, locus_names=target_ids, read_names=ec_ids)

        aln_mat_kallisto.count = counts

        print 'Parsing alignments...'
        widgets = [Bar('>'), ' ', ETA(), ' ', Percentage()]
        pbar = ProgressBar(widgets=widgets, maxval=num_alignments*3).start()
        counter = 0

        for i in xrange(0, num_alignments*3, 3):
            rid = alignments[i]
            lid = alignments[i+1]
            temp_bits = alignments[i+2]

            counter += 1
            pbar.update(i)
            if temp_bits == 0:
                continue

            bits = simple_from_one(temp_bits, num_haplotypes)
            for hid, b in enumerate(bits):
                if b:
                    aln_mat_kallisto.set_value(lid, hid, rid, 1)

        pbar.finish()

        print "Finalizing..."
        aln_mat_kallisto.finalize()
        aln_mat_kallisto.save(emase_file_name, title='KALLISTOALIGN')

        print "DONE"

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", metavar="KallistoAlignFile")
    parser.add_argument("-a", "--apm", dest="apm", metavar="APMFile")

    args = parser.parse_args()

    if not args.input:
        print "No input file specified"
        sys.exit()

    if not args.apm:
        print "No APM file specified"
        sys.exit()

    emasify(args.input, args.apm)

