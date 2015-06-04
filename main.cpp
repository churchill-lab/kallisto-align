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

#include "cxxopts.hpp"

#include <zlib.h>
#include <set>


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
    for (int i = 0; i < b.size(); ++i)
        if (b[i] != 0)
            c |= 1 << i;
    return c;
}

std::vector<int> simple_from_one(int c, int size) {
    std::vector<int> ret;
    for (int i = 0; i < size; ++i)
        ret.push_back((c & (1 << i)) != 0);
    return ret;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

template<typename Index, typename TranscriptCollector>
void createAlignments(Index &index, const ProgramOptions &opt, TranscriptCollector &tc, bool version_ec,
                      std::string emase_binary_file) {

    // need to receive an index map
    std::ios_base::sync_with_stdio(false);

    char delim = '_';

    // read     => idx, read id, sequence
    // target   => idx, target_id
    // hits     => read_idx [target_idx 0,1]

    int tlencount = 10000;
    size_t numreads = 0;
    size_t nummapped = 0;

    // maintarget_idx, haplotype_idx, target_idx
    std::unordered_map<int, std::map<int, int>> map_mainidx_to_all_targets;

    // target_idx, maintarget_idx
    std::unordered_map<int, int> map_targetidx_to_mainidx;

    // maintargetname, maintarget_idx
    std::unordered_map<std::string, int> map_targetname_to_main;

    std::vector<std::string> main_targets;
    std::set<std::string> set_targets;

    std::set<std::string> set_haps;
    std::vector<std::string> haps;
    std::unordered_map<std::string, int> map_haps;


    int idx = 0;
    int main_idx = 0;
    for (auto name : index.target_names_) {
        std::vector<std::string> elems = split(name, delim);

        if (set_targets.insert(elems[0]).second) {
            main_targets.push_back(elems[0]);
            map_targetname_to_main[elems[0]] = main_targets.size() - 1;
        }

        if (set_haps.insert(elems[1]).second) {
            haps.push_back(elems[1]);
            map_haps[elems[1]] = haps.size() - 1;
        }

        main_idx = map_targetname_to_main[elems[0]];
        map_targetidx_to_mainidx[idx] = main_idx;
        map_mainidx_to_all_targets[main_idx][map_haps[elems[1]]] = idx;

        idx++;
    }

    /*

    for (auto s : main_targets) {
        int i = map_targetname_to_main[s];
        std::cout << "TARGET NAME:" << s << ", has main index " << i << std::endl;
        std::map<int, int> c = map_mainidx_to_all_targets[i];

        for (int a = 0; a < haps.size(); a++) {
            std::cout << a << ", "<< haps[a] << ", " << c[a] << std::endl;
        }
    }
     */

    gzFile fp1 = 0, fp2 = 0;
    kseq_t *seq1 = 0, *seq2 = 0;
    std::vector<std::pair<int, int>> v1, v2;
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
    std::cout << PROGNAME << " Pseudo-Aligning..." << std::endl;

    // [readidx] = readid
    std::vector<std::string> read1_ids;
    std::vector<std::string> read2_ids;

    // vector of hash of ints
    std::vector<std::map<int, std::vector<int>>> read1_map_bits;
    std::vector<std::map<int, std::vector<int>>> read2_map_bits;

    // read indices
    int read1_idx = 0;
    int read2_idx = 0;

    std::vector<int> ec_ids;
    // ec_class to index
    std::map<int, int> ec_map;
    std::map<int, std::map<int, std::vector<int>>> ec_bits;

    int ec_index = 0;

    bool ec_write = version_ec;

    int haps_size = haps.size();
    int num_alignment_rows = 0;

    std::vector<int> pos;
    pos.reserve(haps_size);

    std::vector<std::string> elems;
    elems.reserve(2);

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
                    if (ec_map.count(ec) <= 0) {
                        ec_map[ec] = ec_index;
                        ec_ids.push_back(ec);
                        ec_index++;
                    }
                } else {
                    read1_ids.push_back(seq1->name.s);
                    read1_idx++;

                    if (paired) {
                        read2_ids.push_back(seq2->name.s);
                        read2_idx++;
                    }
                }

                // vector elements are index into targets
                const std::vector<int> &vec = index.ecmap[ec];

                //for (auto v : vec) {
                //    std::cout << "v=" << v << ", " << index.target_names_[v] << std::endl;
                //}

                //set_targets.clear();
                std::set<int> int_targets;
                elems.clear();

                // maintarget_idx => bits
                std::map<int, std::vector<int>> _mapper;

                if (ec_write) {
                    if (ec_bits.count(ec) > 0) {
                        _mapper = ec_bits[ec_index - 1];
                    }
                }


                for (int a : vec) {
                    // get the main target
                    int maintarget_idx = map_targetidx_to_mainidx[a];

                    //std::cout << "MAIN TARGET=" << maintarget_idx << ", " << main_targets[maintarget_idx] << std::endl;

                    auto &loc = map_mainidx_to_all_targets[maintarget_idx];

                    if (_mapper.count(maintarget_idx) <= 0) {
                        // maintarget_idx not in mapper
                        for (int h = 0; h < haps_size; h++) {
                            _mapper[maintarget_idx].push_back(0);
                        }
                    }

                    // loop through the haplotypes
                    for (int h = 0; h < haps_size; h++) {
                        //std::cout << "COMPARING " << a << " to " << loc[h] << std::endl;
                        if (a == loc[h]) {
                            _mapper[maintarget_idx][h] = 1;
                            break;
                        }
                    }
                }

                num_alignment_rows += _mapper.size();

                if (ec_write) {
                    ec_bits[ec_index - 1] = _mapper;
                } else {
                    read1_map_bits.push_back(_mapper);
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

    std::cout << PROGNAME << "  Total Targets: " << idx << std::endl;
    std::cout << PROGNAME << "   Main Targets: " << main_targets.size() << std::endl;
    std::cout << PROGNAME << "     Haplotypes: " << haps.size() << std::endl;
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
    int version = 1;

    if (ec_write) {
        std::cout << PROGNAME << " Creating Equivalance Class File..." << std::endl << std::endl;
    } else {
        std::cout << PROGNAME << " Creating Read File..." << std::endl << std::endl;
        version = 0;
    }

    // 1. version
    out.write((char *) &version, sizeof(version));

    // 2. write number of targets
    //std::cout << std::endl << "TARGETS" << std::endl;
    int main_targets_size = main_targets.size();
    out.write((char *) &main_targets_size, sizeof(main_targets_size));

    // 3. write targets
    for (int i = 0; i < main_targets.size(); i++) {
        auto t = main_targets[i];
        tmp_size = strlen(t.c_str());
        out.write((char *) &tmp_size, sizeof(tmp_size));
        out.write(t.c_str(), tmp_size);

        //std::cout << i << ": " << t.c_str() << std::endl;
    }

    // 4. write haplotypes
    //std::cout << std::endl << "HAPLOTYPES" << std::endl;
    int haplotypes_size = haps.size();
    out.write((char *) &haplotypes_size, sizeof(haplotypes_size));

    for (int i = 0; i < haps.size(); i++) {
        auto t = haps[i];
        tmp_size = strlen(t.c_str());
        out.write((char *) &tmp_size, sizeof(tmp_size));
        out.write(t.c_str(), tmp_size);

        //std::cout << i << ": " << t.c_str() << std::endl;
    }

    if (ec_write) {
        // 1. write number of ecids ids
        //std::cout << std::endl << "EC IDS" << std::endl;
        int ec_ids_size = ec_ids.size();
        out.write((char *) &ec_ids_size, sizeof(ec_ids_size));

        for (int id = 0; id < ec_ids.size(); id++) {
            int count = tc.counts[ec_ids[id]];
            out.write((char *) &count, sizeof(count));

            //std::cout << "id=" << id << ", ec = " << ec_ids[id] << ", count = " << count << std::endl;
        }

        // 5. write the ecid, target bits
        int _num_mappings = 0;
        for (int i = 0; i < ec_bits.size(); i++) {
            auto _ecs = ec_bits[i];
            for (auto iterator = _ecs.begin(); iterator != _ecs.end(); iterator++) {
                _num_mappings++;
            }
        }
        out.write((char *) &_num_mappings, sizeof(_num_mappings));

        std::cout << _num_mappings << std::endl;

        for (auto i = ec_bits.begin(); i != ec_bits.end(); i++) {
            for (auto iterator = i->second.begin(); iterator != i->second.end(); iterator++) {
                int idx_ecs = i->first;
                int idx_target = iterator->first;
                int bits = simple_to_one(iterator->second);

                out.write((char *) &idx_ecs, sizeof(idx_ecs));
                out.write((char *) &idx_target, sizeof(idx_target));
                out.write((char *) &bits, sizeof(bits));

                /*
                std::cout << "[" << idx_ecs << "] " << ec_ids[idx_ecs];
                std::cout << "\t"<< tc.counts[ec_ids[idx_ecs]];
                std::cout << "\t[" << idx_target << "] " << main_targets[idx_target];
                std::cout << "\t[" << bits << "] ";

                for (auto b : iterator->second) {
                    std::cout << b;
                }

                std::cout << std::endl;
                */
            }
        }
    } else {
        // 2. write number of read ids
        int reads1_size = read1_ids.size();
        out.write((char *) &reads1_size, sizeof(reads1_size));

        //std::cout << "reads1_size=" << reads1_size << std::endl;

        // 3. write read ids
        for (int i = 0; i < read1_ids.size(); i++) {
            auto r = read1_ids[i];
            tmp_size = strlen(r.c_str());
            out.write((char *) &tmp_size, sizeof(tmp_size));
            out.write(r.c_str(), tmp_size);
        }

        // 7. write read info
        out.write((char *) &num_alignment_rows, sizeof(num_alignment_rows));

        for (int i = 0; i < read1_map_bits.size(); i++) {
            auto rb = read1_map_bits[i];
            for (auto iterator = rb.begin(); iterator != rb.end(); iterator++) {
                int idx_read = i;
                int idx_target = iterator->first;
                int bits = simple_to_one(iterator->second);

                out.write((char *) &idx_read, sizeof(idx_read));
                out.write((char *) &idx_target, sizeof(idx_target));
                out.write((char *) &bits, sizeof(bits));

                /*
                std::cout << "[" << idx_read << "] " << read1_ids[idx_read];
                std::cout << "\t[" << idx_target << "] " << main_targets[idx_target];
                std::cout << "\t[" << bits << "] ";

                for (auto b : iterator->second) {
                    std::cout << b;
                }

                std::cout << std::endl;
                */
            }
        }
    }

    out.flush();
    out.close();
}

/**
 * Useful for looking at the stored binary data.
 */
void loadAlignments(const std::string &emase_binary_file) {
    std::ifstream in;

    in.open(emase_binary_file, std::ios::in | std::ios::binary);

    if (!in.is_open()) {
        std::cerr << PROGNAME << " Error: index input file could not be opened!";
        exit(1);
    }

    std::vector<std::string> reads1;
    std::vector<int> ec_classes_counts;
    std::vector<std::string> haplotypes;
    std::vector<std::string> main_targets;
    std::map<int, std::map<int, std::vector<int>>> reads1_bits;

    int tmp_size;
    int bufsz = 1024;
    char *buffer = new char[bufsz];
    bool ec_read = true;

    // get the version
    int binary_type;
    in.read((char *) &binary_type, sizeof(binary_type));

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

    // 3. read the number of main_targets
    int main_target_size;
    in.read((char *) &main_target_size, sizeof(main_target_size));

    std::cout << PROGNAME << " Number of Main Targets: " << main_target_size << std::endl;

    // 4. read in the main targets
    memset(buffer, 0, bufsz);
    for (auto i = 0; i < main_target_size; ++i) {
        in.read((char *) &tmp_size, sizeof(tmp_size));
        if (tmp_size + 1 > bufsz) {
            delete[] buffer;
            bufsz = 2 * (tmp_size + 1);
            buffer = new char[bufsz];
        }
        memset(buffer, 0, bufsz);
        in.read(buffer, tmp_size);
        main_targets.push_back(std::string(buffer));
    }

    // 5. read in the number of haplotypes
    int hap_size;
    in.read((char *) &hap_size, sizeof(hap_size));

    std::cout << PROGNAME << " Number of Haplotypes: " << hap_size << std::endl;

    // 6. read in the haplotypes
    memset(buffer, 0, bufsz);
    for (auto i = 0; i < hap_size; ++i) {
        in.read((char *) &tmp_size, sizeof(tmp_size));
        if (tmp_size + 1 > bufsz) {
            delete[] buffer;
            bufsz = 2 * (tmp_size + 1);
            buffer = new char[bufsz];
        }
        memset(buffer, 0, bufsz);
        in.read(buffer, tmp_size);
        haplotypes.push_back(std::string(buffer));
    }


    if (ec_read) {
        // 1. read the number of ec classes and counts
        int ec_classes_size;
        in.read((char *) &ec_classes_size, sizeof(ec_classes_size));

        // 2. read in the counts
        for (auto i = 0; i < ec_classes_size; ++i) {
            int idx_count;
            in.read((char *) &idx_count, sizeof(idx_count));
            ec_classes_counts.push_back(idx_count);
        }
        std::cout << PROGNAME << " Number of EC Classes: " << ec_classes_size << std::endl;

        // 7. read the ec class info
        int num_alignments;
        in.read((char *) &num_alignments, sizeof(num_alignments));

        std::cout << PROGNAME << " Number of Alignment Records: " << num_alignments << std::endl;

        for (auto ri = 0; ri < num_alignments; ++ri) {

            int idx_ec;
            in.read((char *) &idx_ec, sizeof(idx_ec));

            int idx_target;
            in.read((char *) &idx_target, sizeof(idx_target));

            int bits;
            in.read((char *) &bits, sizeof(bits));
            std::vector<int> all_bits = simple_from_one(bits, hap_size);

            std::map<int, std::vector<int>> _yn;
            if (reads1_bits.count(idx_ec) > 0) {
                _yn = reads1_bits[idx_ec];
            }

            _yn[idx_target] = all_bits;
            reads1_bits[idx_ec] = _yn;
        }

        for (auto i = reads1_bits.begin(); i != reads1_bits.end(); i++) {
            for (auto j = i->second.begin(); j != i->second.end(); j++) {
                int idx_read = i->first;
                int idx_target = j->first;

                std::cout << idx_read << "\t";
                std::cout << ec_classes_counts[idx_read] << "\t";
                std::cout << main_targets[idx_target] << "\t";

                for (auto &a : j->second) {
                    std::cout << a;
                }

                std::cout << std::endl;
            }
        }

    } else {
        // 1. read the number of reads
        int reads1_size;
        in.read((char *) &reads1_size, sizeof(reads1_size));

        // 2. read in the reads
        for (auto i = 0; i < reads1_size; ++i) {
            in.read((char *) &tmp_size, sizeof(tmp_size));
            if (tmp_size + 1 > bufsz) {
                delete[] buffer;
                bufsz = 2 * (tmp_size + 1);
                buffer = new char[bufsz];
            }
            memset(buffer, 0, bufsz);
            in.read(buffer, tmp_size);
            reads1.push_back(std::string(buffer));
        }
        std::cout << PROGNAME << " Number of reads: " << reads1_size << std::endl;

        // 7. read the read_info
        int num_alignments;
        in.read((char *) &num_alignments, sizeof(num_alignments));

        std::cout << PROGNAME << " Number of Alignment Records: " << num_alignments << std::endl;

        for (auto ri = 0; ri < num_alignments; ++ri) {
            int idx_read;
            in.read((char *) &idx_read, sizeof(idx_read));

            int idx_target;
            in.read((char *) &idx_target, sizeof(idx_target));

            int bits;
            in.read((char *) &bits, sizeof(bits));
            std::vector<int> all_bits = simple_from_one(bits, hap_size);

            std::map<int, std::vector<int>> _yn;
            if (reads1_bits.count(idx_read) > 0) {
                _yn = reads1_bits[idx_read];
            }

            _yn[idx_target] = all_bits;
            reads1_bits[idx_read] = _yn;
        }

        for (auto i = reads1_bits.begin(); i != reads1_bits.end(); i++) {
            for (auto j = i->second.begin(); j != i->second.end(); j++) {
                int idx_read = i->first;
                int idx_target = j->first;

                std::cout << reads1[idx_read] << "\t";
                std::cout << main_targets[idx_target] << "\t";

                for (auto &a : j->second) {
                    std::cout << a;
                }

                std::cout << std::endl;
            }
        }
    }

    delete[] buffer;
    buffer = nullptr;
    in.close();
}


/**
 * Let's do this!
 *
 * V1 = READS
 * V2 = EQUIVALENCE CLASSES (DEFAULT)
 *
 */
int main(int argc, char *argv[]) {
    try {

        bool load = false;
        bool version_read = false;

        // simple argument parsing, nothing special
        cxxopts::Options options(argv[0], " - export kallisto to new binary format");
        options.add_options()
                ("help", "Print help")
                ("l,load", "view the binary file", cxxopts::value<bool>(load))
                ("reads", "create Reads version binary format, else Equivalence Class",
                 cxxopts::value<bool>(version_read))
                ("f,file", "Input Fastq File", cxxopts::value<std::vector<std::string>>(), "FASTQ FILE")
                ("i,index", "Input Index File", cxxopts::value<std::vector<std::string>>(), "INDEX FILE")
                ("b,bin", "Emase Binary File", cxxopts::value<std::string>());

        options.parse(argc, argv);

        if (options.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        // the binary file we are creating or reading
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

            // if not loading, assume creating

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

            // load the index

            KmerIndex index(opt);
            index.load(opt);

            // create the alignments

            MinCollector collection(index, opt);
            createAlignments(index, opt, collection, !version_read, binaryfile);
        }

        std::cout << PROGNAME << " Done" << std::endl << std::endl;

    } catch (const cxxopts::OptionException &e) {
        std::cerr << PROGNAME << " error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    return 0;
}