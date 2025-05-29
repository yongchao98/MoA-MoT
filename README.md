# R1-Code-Interpreter

Our code is based on Llama-factory (https://github.com/hiyouga/LLaMA-Factory) and VeRL (https://github.com/volcengine/verl) for the SFT and RL training and inference on multiple GPUs.
Environment setup:
```
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

GRPO training:
```
In train_grpo_3B.sh, fill your wandb key and python local path in line 1 and line 2.
sh train_grpo_3B.sh
```


```
git clone https://github.com/yongchao98/R1-Code-Interpreter.git
cd R1-Code-Interpreter
conda create -n R1_code_inter python=3.11
conda activate R1_code_inter
pip install reasoning-gym
pip install -r requirements.txt
```

You need to change train_grpo_0_model_base_model.sh so that it requests 8 H100 GPUs or better GPUs. Also modify the python path in the train_grpo_0_model_base_model.sh.
```
cd Search-R1
sh train_grpo_0_model_base_model.sh
```
