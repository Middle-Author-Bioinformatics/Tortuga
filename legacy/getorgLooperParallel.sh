#!/bin/bash

# Set max concurrent jobs (tune based on CPU capacity)
MAX_JOBS=16
job_count=0

for NUM in {0.01,0.05,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90}; do
  (
    SUB=S_${NUM}
    mkdir -p $SUB

    for i in *R1_001.fastq.gz; do
        echo "Processing $i at subsample $NUM"
        r2=${i%_R1_*}_R2_001.fastq.gz

        # Subsample
        seqtk sample $i $NUM > $SUB/${i%.*}
        seqtk sample $r2 $NUM > $SUB/${i%_R1_*}_R2_001.fastq

        # QC and trimming
        ./QC_reads.sh 12 "${i%-SK*}" "$SUB" "$SUB/${i%_R1*}_R1_001.fastq" "$SUB/${i%_R1*}_R2_001.fastq" trimmomatic /home/ark/bin/Trimmomatic-0.39 flash2 pigz 30

        # Organelle assembly
        eval "$(conda shell.bash hook)"
        conda activate getorganelle
        get_organelle_from_reads.py \
            -1 "$SUB/${i%-SK*}.notCombined.final_1P.gz" \
            -2 "$SUB/${i%-SK*}.notCombined.final_2P.gz" \
            -u "$SUB/${i%-SK*}_reads.U.fq.gz" \
            -F animal_mt -o "$SUB/${i%-SKQ23*}-mitos" -t 12 -R 30 \
            --max-reads inf --disentangle-time-limit 200000000 --reduce-reads-for-coverage inf

        # Cleanup
        rm "$SUB/${i%.*}"*
        conda deactivate
    done
  ) &

  ((job_count++))
  # Wait if max parallel jobs are reached
  if [[ $job_count -ge $MAX_JOBS ]]; then
    wait
    job_count=0
  fi
done

# Final wait for all background jobs to complete
wait
