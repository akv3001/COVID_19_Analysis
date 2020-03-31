#!/bin/bash

data_path="/athena/elementolab/scratch/akv3001/COVID_19/COV2_RNASeq_sinai/Raw_Fastq"
data_input="/athena/elementolab/scratch/akv3001/COVID_19/COV2_RNASeq_sinai"

unset http_proxy
unset HTTPS_PROXY
unset FTP_PROXY

for i in $(cat ${data_input}/GSE147507_withGSM | cut -f2|sort |uniq)
        do echo $i

           mkdir ${data_path}/${i}      # Make dir for the main GSMID ( for unique sample)
           SRRall=$(cat ${data_input}/GSE147507_withGSM| grep $i|cut -f1)
           echo $SRRall # grabs all SRR ids's / multiple runs associated with GSM ID

           for eachSRR in $SRRall
                 do echo $eachSRR # Each Run
                    cd  ${data_path}/${i}
                    fastq-dump --split-files -I --gzip $eachSRR
                done
         cd ${data_input}

        done
