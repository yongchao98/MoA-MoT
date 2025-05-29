# R1-Code-Interpreter: Training LLMs to Reason with Code via Supervised and Reinforcement Learning

Our code is based on [Llama-factory](https://github.com/hiyouga/LLaMA-Factory), [VeRL](https://github.com/volcengine/verl), and [Search-R1](https://github.com/PeterGriffinJin/Search-R1?tab=readme-ov-file) for the SFT and RL training and inference on multiple GPUs.

## Environment_Setup
```
git clone https://github.com/yongchao98/R1-Code-Interpreter.git
cd R1-Code-Interpreter
conda create -n R1_code_inter python=3.11
conda activate R1_code_inter
pip install reasoning-gym
pip install -r requirements.txt
pip3 install torch==2.4.0 --index-url https://download.pytorch.org/whl/cu124
pip3 install flash-attn --no-build-isolation
cd verl
pip3 install -e .
cd ../LLaMA-Factory
pip install -e ".[torch,metrics]"
pip install deepspeed==0.15.2
pip install --upgrade huggingface_hub
huggingface-cli login
```

SFT training:
```
cd LLaMA-Factory
sh finetune_qwen_7b_1M.sh
```

GRPO training (In train_grpo_3B.sh, fill your wandb key and python local path in line 1 and line 2):
```
sh train_grpo_3B.sh
```

## Citation
```md
@misc{chen2025r1codeinterpretertrainingllmsreason,
      title={R1-Code-Interpreter: Training LLMs to Reason with Code via Supervised and Reinforcement Learning}, 
      author={Yongchao Chen and Yueying Liu and Junwei Zhou and Yilun Hao and Jingquan Wang and Yang Zhang and Chuchu Fan},
      year={2025},
      eprint={2505.21668},
      archivePrefix={arXiv},
      primaryClass={cs.AI},
      url={https://arxiv.org/abs/2505.21668}, 
}
```
