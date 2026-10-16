#!/bin/bash
# submit 30 chunk jobs (+2°C)

N_CHUNKS=30
WARM=2
SCENARIO=w2

SCRIPT=/storage/brno2/home/urtzien/individuals_DEB/R/3.1_i_europe_chunk.sh

for i in $(seq 1 $N_CHUNKS); do
  JOB="eu_w2_${i}of${N_CHUNKS}"
  echo "Submitting $JOB ..."
  qsub -N $JOB -v CHUNK=$i,N_CHUNKS=$N_CHUNKS,WARM=$WARM,SCENARIO=$SCENARIO $SCRIPT
done

echo "All $N_CHUNKS jobs submitted."
