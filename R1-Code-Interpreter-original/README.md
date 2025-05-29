# R1-Code-Interpreter

python benchmark_test_baseline_new_1.py

In generation_models.py, set the API keys, do not spread to others.

Logic_Game_func.py comprises many environments and one unified function load_task_dataset to load questions and verify_solution_func_gather to check answers.

Task for Yueying: There are some developped tasks not included in this unified Logic_Game_func.py. Tasks: number_multiply, BoxNet1, BoxLift, Blocksworld, big_bench_hard, gsm, math_counting_and_probability, math_geometry. I need you to move these tasks into Logic_Game_func.py. Examples are game24 and letters. Their original environments are in Benchmark dir, run_{env_name}_baseline_methods.py. You need to look into these files and move their codes into Logic_Game_func.py, Note big_bench_hard has four tasks here. So you need to split into 4 task_names.

Task for Junwei: search new tasks in CodeIO and ReasoningGym and then develop tasks as in Logic_Game_func.py. Main components are question loader and answer verifier.
