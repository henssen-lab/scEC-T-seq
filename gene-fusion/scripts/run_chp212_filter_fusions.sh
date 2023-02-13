#!/bin/bash

# filter fusions

# get all breakpoints
cat ../chp212_plate5a6/fusions.tsv | awk '{split($5, a, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"$1}' > ../chp212_plate5a6/fusions_breakpoint_left.txt
cat ../chp212_plate5a6/fusions.tsv | awk '{split($6, a, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"$2}' >> ../chp212_plate5a6/fusions_breakpoint_right.txt

# keep breakpoints on chr2 and sort
tail -n+2 ../chp212_plate5a6/fusions.tsv | awk '{split($5, a, ":");split($6, b, ":"); if(a[1]==2 && b[1] == 2) {print}}' > f.txt
grep -P "^2\t" ../chp212_plate5a6/fusions_breakpoint_left.txt | sort -k1,1 -k2,2n > ../chp212_plate5a6/fusions_breakpoint_left_sorted.txt
grep -P "^2\t" ../chp212_plate5a6/fusions_breakpoint_right.txt | sort -k1,1 -k2,2n > ../chp212_plate5a6/fusions_breakpoint_right_sorted.txt

# overlap with the amplicon regions (+- 10MB)

bedtools intersect -wa -a ../chp212_plate5a6/fusions_breakpoint_left_sorted.txt -b ../chp212_plate5a6/chp212_coordinates_10000000bp.txt | uniq > tmp1
bedtools intersect -wa -a ../chp212_plate5a6/fusions_breakpoint_right_sorted.txt -b ../chp212_plate5a6/chp212_coordinates_10000000bp.txt | uniq > tmp2

head -1 ../chp212_plate5a6/fusions.tsv > ../chp212_plate5a6/fusions_candidates.tsv
cat tmp1 | awk '{print $3}' | grep -f- f.txt > tmp3
cat tmp2 | awk '{print $3}' | grep -f- tmp3 >> ../chp212_plate5a6/fusions_candidates.tsv

rm tmp1 tmp2 tmp3 ../chp212_plate5a6/fusions_breakpoint_left.txt ../chp212_plate5a6/fusions_breakpoint_right.txt ../chp212_plate5a6/fusions_breakpoint_left_sorted.txt ../chp212_plate5a6/fusions_breakpoint_right_sorted.txt

tail -n+2 ../chp212_plate5a6/fusions_candidates.tsv | awk -F"\t" '$13 >= 30 && $14 >= 30 && ($10+$12+$11)>30  &&  2*($10+$12+$11)/($13+$14)>0.2 {print}' > ../chp212_plate5a6/fusions_confident.tsv
tail -n+2 ../chp212_plate5a6/fusions_candidates.tsv | awk -F"\t" '$13 >= 50 && $14 >= 50 &&  2*($10+$12+$11)/($13+$14)>0.3 {print}' >> ../chp212_plate5a6/fusions_confident.tsv
cat ../chp212_plate5a6/fusions_confident.tsv | sort | uniq > tmp

head -1 ../chp212_plate5a6/fusions.tsv > ../chp212_plate5a6/fusions_confident.tsv
cat tmp >> ../chp212_plate5a6/fusions_confident.tsv

tail -n+2 ../chp212_plate5a6/fusions_candidates.tsv | awk '{split($5, a, ":"); split($6, b, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"b[1]"\t"b[2]-1"\t"b[2]}' > ../chp212_plate5a6/fusions_candidates.bedpe
tail -n+2 ../chp212_plate5a6/fusions_confident.tsv | awk '{split($5, a, ":"); split($6, b, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"b[1]"\t"b[2]-1"\t"b[2]}' > ../chp212_plate5a6/fusions_confident.bedpe

head -1 ../chp212_plate5a6/fusions.tsv > ../chp212_plate5a6/fusions_confident_middle.tsv
tail -n+2 ../chp212_plate5a6/fusions_candidates.tsv | awk -F"\t" '$13 >= 30 && $14 >= 30 && ($10+$12+$11)>10 {print}' >> ../chp212_plate5a6/fusions_confident_middle.tsv

# store high
head -1 ../chp212_plate5a6/fusions.tsv > ../chp212_plate5a6/fusions_high.tsv
grep "high" f.txt >> ../chp212_plate5a6/fusions_high.tsv
rm f.txt
