# 200 GB job, 8 GPUs on one host
bsub -J grpo0 \
     -n 4 \
     -gpu "num=8:mode=exclusive_process" \
     -R 'span[hosts=1] rusage[mem=204800]' \
     -M 204800 \
     -o out_R1_2.log \
     -e out_R1_2.err \
     bash train_grpo_0_model.sh

