#!/bin/bash

if [ "$#" -eq 1 ]; then
  export ACCESSION=$1
fi

# Download from SRA
if [ ! -f "${ACCESSION}.download.done" ]; then
  {
    echo "Downloading $ACCESSION"
    fasterq-dump --split-files -p $ACCESSION
  } && {
    touch ${ACCESSION}.download.done
  }
fi

if [ ! -f "${ACCESSION}.download.done" ]; then
  echo "Failed to download."
  exit 1
fi

# Fix the default read names to have /1 and /2 slashes
if [ ! -f "${ACCESSION}.compression.1.done" ]; then
  {
    echo "Compressing/fixing read names ${ACCESSION}_1.fastq"
    sed -r 's/(^[\@\+]SRR\S+)/\1\/1/' ${ACCESSION}_1.fastq | gzip > ${ACCESSION}_1.fastq.gz
    rm ${ACCESSION}_1.fastq
  } && {
    touch ${ACCESSION}.compression.1.done
  }
fi

if [ ! -f "${ACCESSION}.compression.1.done" ]; then
  echo "Failed to compress fastq 1."
  exit 1
fi


if [ ! -f "${ACCESSION}.compression.2.done" ]; then
  {
    echo "Compressing/fixing read names ${ACCESSION}_2.fastq"
    sed -r 's/(^[\@\+]SRR\S+)/\1\/2/' ${ACCESSION}_2.fastq | gzip > ${ACCESSION}_2.fastq.gz
    rm ${ACCESSION}_2.fastq
  } && {
    touch ${ACCESSION}.compression.2.done
  }
fi

if [ ! -f "${ACCESSION}.compression.2.done" ]; then
  echo "Failed to compress fastq 2."
  exit 1
fi

