#!/bin/bash

# get all breakpoints
cat ../tr14_plate3/fusions.tsv | awk '{split($5, a, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"$1}' > ../tr14_plate3/fusions_breakpoint_left.txt
cat ../tr14_plate3/fusions.tsv | awk '{split($6, a, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"$2}' >> ../tr14_plate3/fusions_breakpoint_right.txt

# keep breakpoints on chr2 and sort
tail -n+2 ../tr14_plate3/fusions.tsv | awk '{split($5, a, ":");split($6, b, ":"); if((a[1] == 2 || a[1] == 12) && (b[1] == 2 || b[1] == 12)) {print}}' > f.txt
grep -P "^2\t" ../tr14_plate3/fusions_breakpoint_left.txt | sort -k1,1 -k2,2n > ../tr14_plate3/fusions_breakpoint_left_sorted.txt
grep -P "^12\t" ../tr14_plate3/fusions_breakpoint_left.txt | sort -k1,1 -k2,2n >> ../tr14_plate3/fusions_breakpoint_left_sorted.txt

grep -P "^2\t" ../tr14_plate3/fusions_breakpoint_right.txt | sort -k1,1 -k2,2n > ../tr14_plate3/fusions_breakpoint_right_sorted.txt
grep -P "^12\t" ../tr14_plate3/fusions_breakpoint_right.txt | sort -k1,1 -k2,2n >> ../tr14_plate3/fusions_breakpoint_right_sorted.txt

# overlap with the amplicon regions (+- 10MB)

bedtools intersect -wa -a ../tr14_plate3/fusions_breakpoint_left_sorted.txt -b ../tr14_plate3/tr14_coordinates_10000000bp.txt | uniq > tmp1
bedtools intersect -wa -a ../tr14_plate3/fusions_breakpoint_right_sorted.txt -b ../tr14_plate3/tr14_coordinates_10000000bp.txt | uniq > tmp2

head -1 ../tr14_plate3/fusions.tsv > ../tr14_plate3/fusions_candidates.tsv
cat tmp1 | awk '{print $3}' | grep -f- f.txt > tmp3
cat tmp2 | awk '{print $3}' | grep -f- tmp3 >> ../tr14_plate3/fusions_candidates.tsv


rm tmp1 tmp2 tmp3 ../tr14_plate3/fusions_breakpoint_left.txt ../tr14_plate3/fusions_breakpoint_right.txt ../tr14_plate3/fusions_breakpoint_left_sorted.txt ../tr14_plate3/fusions_breakpoint_right_sorted.txt

tail -n+2 ../tr14_plate3/fusions_candidates.tsv | awk -F"\t" '$13 >= 30 && $14 >= 30 && ($10+$12+$11)>30  &&  2*($10+$12+$11)/($13+$14)>0.2 {print}' > ../tr14_plate3/fusions_confident.tsv
tail -n+2 ../tr14_plate3/fusions_candidates.tsv | awk -F"\t" '$13 >= 50 && $14 >= 50 &&  2*($10+$12+$11)/($13+$14)>0.3 {print}' >> ../tr14_plate3/fusions_confident.tsv
cat ../tr14_plate3/fusions_confident.tsv | sort | uniq > tmp

head -1 ../tr14_plate3/fusions.tsv > ../tr14_plate3/fusions_confident.tsv
cat tmp >> ../tr14_plate3/fusions_confident.tsv

tail -n+2 ../tr14_plate3/fusions_candidates.tsv | awk '{split($5, a, ":"); split($6, b, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"b[1]"\t"b[2]-1"\t"b[2]}' > ../tr14_plate3/fusions_candidates.bedpe
tail -n+2 ../tr14_plate3/fusions_confident.tsv | awk '{split($5, a, ":"); split($6, b, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"b[1]"\t"b[2]-1"\t"b[2]}' > ../tr14_plate3/fusions_confident.bedpe

head -1 ../tr14_plate3/fusions.tsv > ../tr14_plate3/fusions_confident_middle.tsv
tail -n+2 ../tr14_plate3/fusions_candidates.tsv | awk -F"\t" '$13 >= 30 && $14 >= 30 && ($10+$12+$11)>10 {print}' >> ../tr14_plate3/fusions_confident_middle.tsv

head -1 ../tr14_plate3/fusions.tsv > ../tr14_plate3/fusions_high.tsv
grep "high" f.txt >> ../tr14_plate3/fusions_high.tsv
rm f.txt
