#!/bin/bash

export PYTHONPATH="/proj/long-multi/kqian/speech_2":$PYTHONPATH
FORCE_TORCHRUN=1 llamafactory-cli train examples/train_full/qwen_3B_1M_full_sft_ds3_R1_CI-1.yaml
