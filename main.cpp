/*
 * Copyright (c) 2015 The Jackson Laboratory
 *
 * This software was developed by Gary Churchill's Lab at The Jackson
 * Laboratory (see http://research.jax.org/faculty/churchill).
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software. If not, see <http://www.gnu.org/licenses/>.
 */

#include "kseq.h"
#include "KmerIndex.h"
#include "MinCollector.h"
#include "common.h"
#include "Kmer.hpp"

#include "cxxopts.hpp"

#include <zlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>


#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

#define PROGNAME "[kallisto-export]"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

int simple_to_one(const std::vector<int> b) {
    int c = 0;
    for (int i=0; i < b.size(); ++i)
        if (b[i] != 0)
            c |= 1 << i;
    return c;
}

std::vector<int> simple_from_one(int c, int size) {
    std::vector<int> ret;
    for (int i=0; i < size; ++i)
        ret.push_back((c & (1<<i)) != 0);
    return ret;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

struct vec_ints {
    std::vector<int> bits;
};
typedef std::map<int, vec_ints> map_ints;

template<typename Index, typename TranscriptCollector>
void createAlignments(Index& index, const ProgramOptions& opt, TranscriptCollector& tc, bool version_ec, std::string emase_binary_file) {

    // need to receive an index map
    std::ios_base::sync_with_stdio(false);

    char delim = '_';

    // read     => idx, read id, sequence
    // target   => idx, target_id
    // hits     => read_idx [target_idx 0,1]

    int tlencount = 10000;
    size_t numreads = 0;
    size_t nummapped = 0;

    // create the lookup for targets
    // simple_id, haplotype, target_idx
    std::map <std::string, std::map<std::string, int>> lookup_targets;

    std::set <std::string> set_haps;
    std::set <std::string> set_targets;
    std::map <std::string, int> main_targets_rev;
    std::vector <std::string> main_targets;
    std::vector <std::string> haps;

    int idx = 0;
    for (auto name : index.target_names_) {
        std::vector <std::string> elems = split(name, delim);

        if (set_haps.count(elems[1]) <= 0) {
            set_haps.insert(elems[1]);
            haps.push_back(elems[1]);
        }

        if (set_targets.count(elems[0]) <= 0) {
            set_targets.insert(elems[0]);
            main_targets.push_back(elems[0]);
            main_targets_rev[elems[0]] = main_targets.size() - 1;
        }

        lookup_targets[elems[0]][elems[1]] = idx++;
    }


    gzFile fp1 = 0, fp2 = 0;
    kseq_t *seq1 = 0, *seq2 = 0;
    std::vector <std::pair<int, int>> v1, v2;
    v1.reserve(1000);
    v2.reserve(1000);

    bool paired = !opt.single_end;

    int l1 = 0, l2 = 0; // length of read

    if (paired) {
        std::cout << PROGNAME << " Mode: paired-end" << std::endl;
    } else {
        std::cout << PROGNAME << " Mode: single-end" << std::endl;
    }

    for (int i = 0; i < opt.files.size(); i += (paired) ? 2 : 1) {
        if (paired) {
            std::cout << PROGNAME << " Processing pair " << (i / 2 + 1) << ": " << opt.files[i] << std::endl
            << "                        " << opt.files[i + 1] << std::endl;
        } else {
            std::cout << PROGNAME << " Processing " << i + 1 << ": " << opt.files[i] << std::endl;
        }
    }

    // for each file
    std::cout << PROGNAME << " Aligning..." << std::endl;

    // [readidx] = <readid, sequence>
    std::vector <std::pair<std::string, std::string>> reads1;
    std::vector <std::pair<std::string, std::string>> reads2;

    // vector of hash of ints
    std::map<int, map_ints> reads1_bits;
    std::map<int, map_ints> reads2_bits;
    std::set <std::string> _reads1;
    std::set <std::string> _reads2;
    int read1_idx = 0;
    int read2_idx = 0;
    int num_alignment_rows = 0;

    std::map<int, map_ints> ec_bits;
    std::set<int> _ecs;
    int ec_index = 0;
    std::vector <int> ec_ids;
    bool ec_write = version_ec;

    for (int i = 0; i < opt.files.size(); i += (paired) ? 2 : 1) {

        fp1 = gzopen(opt.files[i].c_str(), "r");
        seq1 = kseq_init(fp1);
        if (paired) {
            fp2 = gzopen(opt.files[i + 1].c_str(), "r");
            seq2 = kseq_init(fp2);
        }

        // for each read
        while (true) {
            l1 = kseq_read(seq1);
            if (paired) {
                l2 = kseq_read(seq2);
            }
            if (l1 <= 0) {
                break;
            }
            if (paired && l2 <= 0) {
                break;
            }

            numreads++;
            v1.clear();
            v2.clear();
            // process read
            index.match(seq1->seq.s, seq1->seq.l, v1);
            if (paired) {
                index.match(seq2->seq.s, seq2->seq.l, v2);
            }

            // collect the target information
            int ec = tc.collect(v1, v2, !paired);
            if (ec != -1) {
                nummapped++;

                if (ec_write) {
                    if (_ecs.count(ec) <= 0) {
                        _ecs.insert(ec);
                        ec_index++;
                        ec_ids.push_back(ec);
                    }
                } else {
                    if (_reads1.count(seq1->name.s) <= 0) {
                        //std::cout << "R_E_A_D adding : " << read1_idx << " : " << seq1->name.s << std::endl;
                        _reads1.insert(seq1->name.s);
                        reads1.push_back({seq1->name.s, seq1->seq.s});
                        read1_idx++;
                    }

                    if (paired) {
                        if (_reads2.count(seq2->name.s) <= 0) {
                            _reads2.insert(seq2->name.s);
                            reads2.push_back({seq2->name.s, seq2->seq.s});
                            read2_idx++;
                        }
                    }
                }

                // vector elements are index into target names
                const std::vector <int> &vec = index.ecmap[ec];

                //std::vector< std::map<int, std::vector<int>> target_idxs;
                set_targets.clear();

                // get the main transcript id so we can see if all haplotypes match

                map_ints _mapper;

                for (int a : vec) {
                    std::vector <std::string> elems = split(index.target_names_[a], delim);
                    if (set_targets.count(elems[0]) <= 0) {
                        set_targets.insert(elems[0]);

                        std::vector <int> pos;

                        for (auto h : haps) {
                            int target_index = lookup_targets[elems[0]][h];
                            pos.push_back(target_index);
                        }
                        _mapper[main_targets_rev[elems[0]]].bits = pos;
                    }

                }
                num_alignment_rows += set_targets.size();

                map_ints _yn;
                if (reads1_bits.count(read1_idx-1) > 0) {
                    _yn = reads1_bits[read1_idx-1];
                }
                for (auto iterator = _mapper.begin(); iterator != _mapper.end(); iterator++) {
                    std::vector <int> pos;

                    for (auto a : iterator->second.bits) {
                        bool found = false;

                        for (int b : vec) {
                            if (a == b) {
                                found = true;
                                break;
                            }
                        }
                        pos.push_back(found ? 1 : 0);
                    }
                    _yn[iterator->first].bits = pos;
                }

                if (ec_write) {
                    ec_bits[ec_index - 1] = _yn;
                } else {
                    reads1_bits[read1_idx - 1] = _yn;
                }
            }

            if (paired && 0 <= ec && ec < index.num_trans && tlencount > 0) {
                bool allSame = (v1[0].first == ec && v2[0].first == ec) && (v1[0].second == 0 && v2[0].second == 0);

                if (allSame) {
                    // try to map the reads
                    int tl = index.mapPair(seq1->seq.s, seq1->seq.l, seq2->seq.s, seq2->seq.l, ec);
                    if (0 < tl && tl < tc.flens.size()) {
                        tc.flens[tl]++;
                        tlencount--;
                    }
                }
            }

        }
        gzclose(fp1);
        if (paired) {
            gzclose(fp2);
        }
    }

    kseq_destroy(seq1);
    if (paired) {
        kseq_destroy(seq2);
    }

    std::cout << PROGNAME << "     Haplotypes: " << haps.size() << std::endl;
    std::cout << PROGNAME << "  Total Targets: " << idx << std::endl;
    std::cout << PROGNAME << "   Main Targets: " << main_targets.size() << std::endl;
    std::cout << PROGNAME << "          Reads: " << numreads << std::endl;
    std::cout << PROGNAME << "   Mapped Reads: " << nummapped << std::endl;
    std::cout << PROGNAME << "     Alignments: " << num_alignment_rows << std::endl;

    std::ofstream out;
    out.open(emase_binary_file, std::ios::out | std::ios::binary);

    if (!out.is_open()) {
        std::cerr << PROGNAME << " Error: align output file could not be opened!";
        exit(1);
    }

    int tmp_size = 0;

    if (ec_write) {
        std::cout << PROGNAME << " Creating Equivalance Class File..." << std::endl << std::endl;

        // 1. write version
        int version = 1;
        out.write((char *) &version, sizeof(version));


        // 1. write number of ecids ids
        int ec_ids_size = ec_ids.size();
        out.write((char *) &ec_ids_size, sizeof(ec_ids_size));

        for (int id = 0; id < ec_ids.size(); id++) {
            int count = tc.counts[ec_ids[id]];
            out.write((char *) &count, sizeof(count));
        }

        // 2. write number of targets
        int main_targets_size = main_targets.size();
        out.write((char *) &main_targets_size, sizeof(main_targets_size));

        // 3. write targets
        for (int i = 0; i < main_targets.size(); i++) {
            auto t = main_targets[i];
            tmp_size = strlen(t.c_str());
            out.write((char *) &tmp_size, sizeof(tmp_size));
            out.write(t.c_str(), tmp_size);
        }

        // 4. write haplotypes
        int haplotypes_size = haps.size();
        out.write((char *) &haplotypes_size, sizeof(haplotypes_size));

        for (int i = 0; i < haps.size(); i++) {
            auto t = haps[i];
            tmp_size = strlen(t.c_str());
            out.write((char *) &tmp_size, sizeof(tmp_size));
            out.write(t.c_str(), tmp_size);
        }

        // 5. write the ecid, target bits
        int _num_mappings = 0;
        for (auto i = ec_bits.begin(); i != ec_bits.end(); i++) {
            for (auto iterator = i->second.begin(); iterator != i->second.end(); iterator++) {
                _num_mappings++;
            }
        }
        out.write((char *) &_num_mappings, sizeof(_num_mappings));

        for (auto i = ec_bits.begin(); i != ec_bits.end(); i++) {
            for (auto iterator = i->second.begin(); iterator != i->second.end(); iterator++) {
                int idx_ecs = i->first;
                int idx_target = iterator->first;
                int bits = simple_to_one(iterator->second.bits);

                out.write((char *) &idx_ecs, sizeof(idx_ecs));
                out.write((char *) &idx_target, sizeof(idx_target));
                out.write((char *) &bits, sizeof(bits));
            }
        }

    } else {
        std::cout << PROGNAME << " Creating Read File..." << std::endl << std::endl;

        // 1. write version
        int version = 0;
        out.write((char *) &version, sizeof(version));

        // 2. write number of read ids
        int reads1_size = reads1.size();
        out.write((char *) &reads1_size, sizeof(reads1_size));

        // 3. write read ids
        for (int i = 0; i < reads1.size(); i++) {
            auto r = reads1[i];
            tmp_size = strlen(r.first.c_str());
            out.write((char *) &tmp_size, sizeof(tmp_size));

            out.write(r.first.c_str(), tmp_size);
        }

        // 4. write number of targets
        int main_targets_size = main_targets.size();
        out.write((char *) &main_targets_size, sizeof(main_targets_size));

        // 5. write targets
        for (int i = 0; i < main_targets.size(); i++) {
            auto t = main_targets[i];
            tmp_size = strlen(t.c_str());
            out.write((char *) &tmp_size, sizeof(tmp_size));
            out.write(t.c_str(), tmp_size);
        }

        // 6. write haplotypes
        int haplotypes_size = haps.size();
        out.write((char *) &haplotypes_size, sizeof(haplotypes_size));

        for (int i = 0; i < haps.size(); i++) {
            auto t = haps[i];
            tmp_size = strlen(t.c_str());
            out.write((char *) &tmp_size, sizeof(tmp_size));
            out.write(t.c_str(), tmp_size);
        }

        // 7. write read info
        out.write((char *) &num_alignment_rows, sizeof(num_alignment_rows));

        for (auto i = reads1_bits.begin(); i != reads1_bits.end(); i++) {
            for (auto iterator = i->second.begin(); iterator != i->second.end(); iterator++) {
                int idx_read = i->first;
                int idx_target = iterator->first;
                int bits = simple_to_one(iterator->second.bits);

                out.write((char *) &idx_read, sizeof(idx_read));
                out.write((char *) &idx_target, sizeof(idx_target));
                out.write((char *) &bits, sizeof(bits));
            }
        }
    }

    out.flush();
    out.close();

}

void loadAlignments(const std::string& emase_binary_file) {
    std::ifstream in;

    in.open(emase_binary_file, std::ios::in | std::ios::binary);

    if (!in.is_open()) {
        // TODO: better handling
        std::cerr << PROGNAME << " Error: index input file could not be opened!";
        exit(1);
    }

    std::vector<std::string> reads1;
    std::vector<int> ec_classes_counts;
    std::vector<std::string> haplotypes;
    std::vector<std::string> main_targets;
    std::map<int, map_ints > reads1_bits;

    int tmp_size;
    int bufsz = 1024;
    char *buffer = new char[bufsz];
    bool ec_read = true;

    int binary_type;
    in.read((char *)&binary_type, sizeof(binary_type));

    if (binary_type == 1) {
        ec_read = true;
        std::cout << PROGNAME << " File contains Equivalent Classes" << std::endl;
    } else if (binary_type == 0) {
        ec_read = false;
        std::cout << PROGNAME << " File contains Reads" << std::endl;
    } else {
        std::cerr << PROGNAME << " Error: unknown file type or version";
        exit(1);
    }

    if (ec_read) {
        // 1. read the number of ec classes and counts
        int ec_classes_size;
        in.read((char *)&ec_classes_size, sizeof(ec_classes_size));

        // 2. read in the reads
        for (auto i = 0; i < ec_classes_size; ++i) {
            int idx_count;
            in.read((char *) &idx_count, sizeof(idx_count));
            ec_classes_counts.push_back(idx_count);
        }
        std::cout << PROGNAME << " Number of EC Classes: " << ec_classes_size << std::endl;
    } else {
        // 1. read the number of reads
        int reads1_size;
        in.read((char *)&reads1_size, sizeof(reads1_size));

        // 2. read in the reads
        for (auto i = 0; i < reads1_size; ++i) {
            in.read((char *)&tmp_size, sizeof(tmp_size));
            if (tmp_size +1 > bufsz) {
                delete[] buffer;
                bufsz = 2*(tmp_size+1);
                buffer = new char[bufsz];
            }
            memset(buffer,0,bufsz);
            in.read(buffer, tmp_size);
            reads1.push_back(std::string( buffer ));
        }
        std::cout << PROGNAME << " Number of reads: " << reads1_size << std::endl;
    }


    // 3. read the number of main_targets
    int main_target_size;
    in.read((char *)&main_target_size, sizeof(main_target_size));

    std::cout << PROGNAME << " Number of Main Targets: " << main_target_size << std::endl;

    // 4. read in the main targets
    memset(buffer,0,bufsz);
    // 2. read in the reads
    for (auto i = 0; i < main_target_size; ++i) {
        in.read((char *)&tmp_size, sizeof(tmp_size));
        if (tmp_size +1 > bufsz) {
            delete[] buffer;
            bufsz = 2*(tmp_size+1);
            buffer = new char[bufsz];
        }
        memset(buffer,0,bufsz);
        in.read(buffer, tmp_size);
        main_targets.push_back(std::string( buffer ));
    }

    // haplotypes
    int hap_size;
    in.read((char *)&hap_size, sizeof(hap_size));

    std::cout << PROGNAME << " Number of Haplotypes: " << hap_size << std::endl;

    memset(buffer,0,bufsz);
    for (auto i = 0; i < hap_size; ++i) {
        in.read((char *)&tmp_size, sizeof(tmp_size));
        if (tmp_size +1 > bufsz) {
            delete[] buffer;
            bufsz = 2*(tmp_size+1);
            buffer = new char[bufsz];
        }
        memset(buffer,0,bufsz);
        in.read(buffer, tmp_size);
        haplotypes.push_back(std::string( buffer ));
    }

    if (ec_read) {

        // 5. read the ec class info
        int num_alignments;
        in.read((char *)&num_alignments, sizeof(num_alignments));

        std::cout << PROGNAME << " Number of Alignment Records: " << num_alignments << std::endl;

        for (auto ri = 0; ri < num_alignments; ++ri) {

            int idx_ec;
            in.read((char *) &idx_ec, sizeof(idx_ec));

            int idx_target;
            in.read((char *) &idx_target, sizeof(idx_target));

            int bits;
            in.read((char *) &bits, sizeof(bits));
            std::vector<int> all_bits = simple_from_one(bits, hap_size);

            map_ints _yn;
            if (reads1_bits.count(idx_ec) > 0) {
                _yn = reads1_bits[idx_ec];
            }

            _yn[idx_target].bits = all_bits;
            reads1_bits[idx_ec] = _yn;
        }

        for(auto i = reads1_bits.begin(); i != reads1_bits.end(); i++) {
            for(auto j = i->second.begin(); j != i->second.end(); j++) {
                int idx_read = i->first;
                int idx_target = j->first;

                std::cout << idx_read << "\t";
                std::cout << ec_classes_counts[idx_read] << "\t";
                std::cout << main_targets[idx_target] << "\t";

                for (auto &a : j->second.bits) {
                    std::cout << a;
                }

                std::cout << std::endl;
            }
        }
    } else {
        // 5. read the read_info
        int num_alignments;
        in.read((char *)&num_alignments, sizeof(num_alignments));

        std::cout << PROGNAME << " Number of Alignment Records: " << num_alignments << std::endl;

        for (auto ri = 0; ri < num_alignments; ++ri) {

            int idx_read;
            in.read((char *) &idx_read, sizeof(idx_read));

            int idx_target;
            in.read((char *) &idx_target, sizeof(idx_target));

            int bits;
            in.read((char *) &bits, sizeof(bits));
            std::vector<int> all_bits = simple_from_one(bits, hap_size);

            map_ints _yn;
            if (reads1_bits.count(idx_read) > 0) {
                _yn = reads1_bits[idx_read];
            }

            _yn[idx_target].bits = all_bits;
            reads1_bits[idx_read] = _yn;
        }

        for(auto i = reads1_bits.begin(); i != reads1_bits.end(); i++) {
            for(auto j = i->second.begin(); j != i->second.end(); j++) {
                int idx_read = i->first;
                int idx_target = j->first;

                std::cout << reads1[idx_read] << "\t";
                std::cout << main_targets[idx_target] << "\t";

                for (auto &a : j->second.bits) {
                    std::cout << a;
                }

                std::cout << std::endl;
            }
        }
    }

    // delete the buffer
    delete[] buffer;
    buffer=nullptr;

    in.close();
}


/**
 * Let's do this!
 *
 * V1 = READS
 * V2 = EQUIVALENCE CLASSES (DEFAULT)
 *
 */
int main(int argc, char* argv[]) {
    try  {
        cxxopts::Options options(argv[0], " - export kallisto to new binary format");

        bool load = false;
        bool version_read = false;

        options.add_options()
                ("help", "Print help")
                ("l,load", "view the binary file", cxxopts::value<bool>(load))
                ("reads", "create Reads version binary format, else Equivalence Class", cxxopts::value<bool>(version_read))
                ("f,file", "Input Fastq File", cxxopts::value<std::vector<std::string>>(), "FASTQ FILE")
                ("i,index", "Input Index File", cxxopts::value<std::vector<std::string>>(), "INDEX FILE")
                ("b,bin", "Emase Binary File", cxxopts::value<std::string>());

        options.parse(argc, argv);

        if (options.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        std::string binaryfile;

        if (options.count("b")) {
            binaryfile = options["b"].as<std::string>();
        } else {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        if (load) {

            std::cout << PROGNAME << " Loading " << binaryfile << "..." << std::endl;
            loadAlignments(binaryfile);

        } else {

            ProgramOptions opt;
            if (options.count("f")) {
                opt.files = options["f"].as<std::vector<std::string>>();
            } else {
                std::cout << "-f must be supplied" << std::endl;
                std::cout << options.help() << std::endl;
                exit(0);
            }

            if (options.count("i")) {
                opt.index = options["i"].as<std::vector<std::string>>()[0];
            } else {
                std::cout << "-i must be supplied" << std::endl;
                std::cout << options.help() << std::endl;
                exit(0);
            }

            opt.single_end = true;

            std::cout << PROGNAME << " Creating " << binaryfile << "..." << std::endl;
            KmerIndex index(opt);
            index.load(opt);
            MinCollector collection(index, opt);
            createAlignments(index, opt, collection, !version_read, binaryfile);
        }

        std::cout << PROGNAME << " Done" << std::endl << std::endl;

    } catch (const cxxopts::OptionException& e) {
        std::cerr << PROGNAME << " error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    return 0;
}