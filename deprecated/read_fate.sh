#!/bin/bash

set -x

cd /home/groups/hpcbio/projects/LifeSpan/exome-March2016/src/logs

truncate -s -0 ReadFateStats.txt

echo "Lets extract stats from the alignment logs"

grep "tot_reads=" *align*.in > $ReadFateStats.txt
grep "tot_mapped=" *align*.in >> $ReadFateStats.txt
grep "tot_dups=" *align*.in >> $ReadFateStats.txt
grep "perc_dup=" *align*.in >> $ReadFateStats.txt
grep "perc_mapped=" *align*.in >> $ReadFateStats.txt


echo "Lets extract stats from the merge logs"

grep "alignments merged." *merge*.in >> $ReadFateStats.txt

