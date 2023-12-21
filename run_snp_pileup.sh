#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-6:00
#SBATCH -p short
#SBATCH --mem 4000
#SBATCH -o snp_pileup.out
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

# install the snp-pileup software on o2 in a conda environment called "snp_pileup" and then activate the environment in this script
# you may have to edit your .bashrc/.bash_profile file so that the miniconda module is in your $PATH every time you logon to o2.
conda activate snp_pileup

if test -f OV1_Normal2_ligated_cleaned.pileup; then
    echo "File exists."
else
    echo "Does not exist. Running snp-pileup."
    snp-pileup glimpse/OV1_Normal2_ligated_cleaned.bcf OV1_Normal2_ligated_cleaned.pileup -q 10 -Q 20 -P 0 \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Ap1a_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Ap1b_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Col1_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Col2_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Col5-A_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Col6-A_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Dia1-A_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Ile1_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Ile2_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Ld2-A_recal.bam \
        /home/alg2264/data/alex/azenta/30-922534314_OV1_2_normals_scraped/preprocessing/Normal2_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Ov1_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/Per1c_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/PT1_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/PT2_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/PT3_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/PT4_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/PT5_recal.bam \
        /home/alg2264/data/alex/azenta/30-859115236_OV1_ALK1_no_extraction/OV1/preprocessing/PT6_recal.bam
fi

gzip OV1_Normal2_ligated_cleaned.pileup
