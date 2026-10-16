#!/bin/bash
#PBS -N ${JOB_NAME:-i_eu_chunk}
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=10:mem=50gb:scratch_local=50gb
#PBS -m abe
#PBS -M uenriquez@proton.me
#PBS -v CHUNK,N_CHUNKS,WARM,SCENARIO

export TMPDIR=$SCRATCHDIR

# repo on storage
REPO=/storage/brno2/home/urtzien/individuals_DEB

# output
MYOUT=$REPO/results
mkdir -p $MYOUT

# modules
module add r/
module add udunits/
module add proj/
module add geos/
module add gdal/
module add openssl/
export R_LIBS_USER=/storage/brno2/home/urtzien/Rpackages

# copy repo to scratch
cd $SCRATCHDIR || exit 2
cp -r $REPO/* .

# run
cd R/
Rscript 3.1_run_europe_chunk.R --chunk=${CHUNK}/${N_CHUNKS} --warm=${WARM} --scenario=${SCENARIO} \
  > europe_chunk_${SCENARIO}_${CHUNK}of${N_CHUNKS}.log 2>&1

# copy results back
OUTFILE=europe_results_${SCENARIO}_chunk${CHUNK}of${N_CHUNKS}.rds
cp europe_chunk_${SCENARIO}_${CHUNK}of${N_CHUNKS}.log ../data/${OUTFILE} $MYOUT/
