#!/usr/bin/env python
import os.path
import argparse
import re
from statistics import mean
from typing import List

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=" Given a path with CSV files, calculate average aligned reads, total uniquely mapped and HTSeq features")

    parser.add_argument( '-c', '--csv_files',
                            dest= "csv_files",
                            action= "store",
                            default= os.getcwd(),
                            help="Input directory with CSV files from HTSeq")

    parser.add_argument( '-s', '--star_files',
                            dest= "star_files",
                            action= "store",
                            default= os.getcwd(),
                            help="Input directory with log files from STAR ")
    
    parser.add_argument( '-o', '--outpit_dir',
                            dest= "output_dir",
                            action= "store",
                            default= os.getcwd(),
                            help="Input directory to store the results ")


    parser.add_argument('-v', '--verbose',
                        dest="verbose",
                        action="store_true",
                        default=False,
                        help="Print progression log to standard error")

    options = parser.parse_args()

    # CAPTURING THE INPUT FILE(s)

    csv_path = options.csv_files
    star_path = options.star_files
    output_dir = options.output_dir

    csv_listfiles = [ os.path.join(csv_path, f) for f in os.listdir(csv_path) if f.endswith(".csv")]

    star_listfiles = [ os.path.join(star_path, f) for f in os.listdir(star_path) if f.endswith(".final.out")]

    htseq = {}
    star = {}

    for csv_file in csv_listfiles:
        fd = open(csv_file, "r")

        htseq[csv_file] = {}
        htseq[csv_file]["feature"] = 0
        for line in fd:
            if "__" in line:
                split_line = line.strip().split("\t")
                htseq[csv_file][split_line[0][2:]] = int(split_line[1])
            else:
                split_line = line.strip().split("\t")
                htseq[csv_file]["feature"] += int(split_line[1])
        fd.close()

    for star_file in star_listfiles:
        fd = open(star_file, "r")

        star[star_file] = []
        for line in fd:
            if "Uniquely mapped reads number" in line:
                sp_line = line.strip().split("\t")
                star[star_file].append(sp_line[1])
            elif "Uniquely mapped reads %" in line:
                line3 = line.strip().split("\t")
                star[star_file].append(line3[1].strip("%"))
        fd.close()

    samples = []

    for file in htseq.keys():
        filename = re.match(r'.*/(.*)_(L00.*)_htseq.csv', file)
        samples.append(filename.group(1))

    final_star_dict = {}
    final_htseq_dict = {}
    htseq_values = {}

    samples = set(samples)
    for sample in samples:
        final_star_dict[sample] = []
        final_htseq_dict[sample] = {}

    for file, value in htseq.items():
        for sample in samples:
            for sub, number in value.items():
                final_htseq_dict[sample][sub] = []

    for file, values in star.items():
         for sample in samples:
             if sample in file:
                 final_star_dict[sample].append(tuple(values))

    for file, values in htseq.items():
         for sample in samples:
             if sample in file:
                 for key, number in values.items():
                    final_htseq_dict[sample][key].append(number)

    for key, value in final_star_dict.items():
        total_uniq_mapped = []  # type: List[int]
        average_mapped = []
        for comb in value:
            total_uniq_mapped.append(int(comb[0]))
            average_mapped.append(float(comb[1]))
        final_star_dict[key] = [sum(total_uniq_mapped), mean(average_mapped)]

    for key, value in final_star_dict.items():
        print(key, value)

    output_file_content = ""
    htseq_names = ["feature", "ambiguous", "no_feature", "too_low_aQual", "not_aligned"]

    for key, value in final_htseq_dict.items():
        htseq_values[key] = ""
        total_reads = 0
        for sub, number in value.items():
            final_htseq_dict[key][sub] = sum(number)
            total_reads += sum(number)
        print(key, value)
        for key2, value2 in final_star_dict.items():
            if key == key2:
                for name in htseq_names:
                    htseq_values[key] += "%s:%.2f," %(name, value[name]/((value2[0]))*100)
                output_file_content += "%10s\t%11.3f\t%18d\t%5s\n" %(key, value2[1], value2[0], htseq_values[key].strip(","))


    ofd = open("%s/summary_star_htseq.txt" %(output_dir), "w")


    ofd.write("%10s\t%13s\t%20s\t%20s\n" %("Sample", "% Alignment", "Uniquely mapped", "HTSeq summary (%)"))
    ofd.write(output_file_content)
    ofd.close()

