#!/bin/bash

export PYTHONPATH="/proj/long-multi/kqian/speech_2":$PYTHONPATH
FORCE_TORCHRUN=1 llamafactory-cli train examples/train_full/llama3_8B_full_sft_ds3_R1_CI-2.yaml
