# R1-Code-Interpreter

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
