import json
import re
import pandas as pd
import os
import subprocess
import sys
from openai import OpenAI
from generation_models import message_construct_func, message_construct_llama_func, GPT_response, count_total_tokens, extract_code, extract_and_check, LLM_answer_code_checker, save_file_func, paraphrase_with_GPT4
import random
import math
import json
from typing import List, Tuple, Dict
import time
import numpy as np
from prompt import *
from argparse import ArgumentParser

from benchmark.run_Logic_Game_baseline_methods_new import run_logic_game_baselines
'''
from benchmark.run_game24_baseline_methods import run_game24_baselines
from benchmark.run_path_plan_continuous_baseline_methods import run_path_plan_baselines
from benchmark.run_letters_baseline_methods import run_letters_baselines
from benchmark.run_number_multiply_baseline_methods import run_number_multiply_baselines
from benchmark.run_BoxNet1_baseline_methods import run_boxnet1_baselines
from benchmark.run_BoxLift_baseline_methods import run_boxlift_baselines
from benchmark.run_Blocksworld_baseline_methods import run_blocksworld_baselines
from benchmark.run_big_bench_hard_baseline_methods import run_big_bench_hard_baselines
from benchmark.run_gsm_baseline_methods import run_gsm_baselines
from benchmark.run_math_counting_and_probability_baseline_methods import run_MATH_counting_and_probability_baselines
from benchmark.run_math_geometry_baseline_methods import run_MATH_geometry_baselines
from benchmark.run_Logic_Game_CodeSteer_DPO_data_creation import run_logic_game_DPO_data
'''

if __name__ == '__main__':
    # gpt-4o, gpt-4o-mini, gpt-3.5-turbo for OpenAi API

    def log_run_info(log_file, run_info):
        with open(log_file, 'a') as f:
            f.write(run_info + "\n")

    model_name = 'gpt-4o'  # Qwen2.5-14B-1M, Qwen2.5-7B-1M, DeepSeek-R1-8B
    # o3-mini-2025-01-31, gpt-4o, gpt-3.5-turbo, "claude-3-sonnet-20240229", o1, o1-preview, gpt-4o, qwen2_QwQ_32B-preview, 'DeepSeek-R1'

    #args_path = '/n/vlassak_lab/Lab/simulation_data/NLP_robotics/experiment/T5/large_model/llama3/LLaMA-Factory/examples/inference/qwen2_7B.yaml'
    #args_path = '/n/vlassak_lab/Lab/simulation_data/NLP_robotics/experiment/T5/large_model/llama3/LLaMA-Factory/examples/inference/deepseek_r1_8B.yaml'
    args_path = '/n/vlassak_lab/Lab/simulation_data/NLP_robotics/experiment/T5/large_model/llama3/LLaMA-Factory/examples/inference/qwen2_14B.yaml'

    gather_save_input_dir = 'results_gather'
    baseline_method_name = 'All_code_CoT' # 1_only_ques, code_interpreter, AutoGen, All_code_CoT, All_text
    max_sample_num = 5
    ### Test baseline methods
    env_name_list1 = ['eight_queens', 'pooling', 'reversi', 'light_puzzles', 'new_operator', 'mahjong_pattern']

    env_name_list2 = ['statistical_counting', 'synthesis_decomposition', '2048', 'matrix_transformation', 'pattern_recognition', 'constrained_linear_arrangement']

    env_name_list3 = ['string_synthesis', 'logic_puzzle', 'string_insertion', 'letter_logic_diagram', 'standard_sudoku', 'string_deletion_and_modification']

    env_name_list4 = ['string_splitting', 'cryptanalysis', 'permutations_and_combinations', 'logical_equation', 'combinatorial_calculation']

    env_name_list5 = ['game24', 'letters']

    base_path = 'results_gather'
    runtime_list = []
    for task_name in ['letters']:
        start_time = time.time()
        run_logic_game_baselines(task_name, gather_save_input_dir, model_name, baseline_method_name, args_path, max_sample_num)
        end_time = time.time()
        runtime = end_time - start_time
        runtime_list.append(runtime/max_sample_num)

        output_path = base_path + f'/Cost_runtime_gather_{model_name}_{baseline_method_name}.txt'
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'a') as f:
            f.write(f"{runtime/max_sample_num}\n")
        print(f'Mean time cost: {np.mean(runtime_list)}')
        print(f'\nDataset saved to: {output_path}')

    '''
    save_input_dir = 'results_gather/MATH_geometry'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_MATH_geometry_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)

    save_input_dir = 'results_gather/MATH_c_p'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_MATH_counting_and_probability_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)

    save_input_dir = 'results_gather/gsm'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_gsm_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)

    save_input_dir = 'results_gather/BIG-Bench-Hard'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_big_bench_hard_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)
    '''

    '''
    save_input_dir = 'results_gather/BoxNet1'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_boxnet1_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)
    
    save_input_dir = 'results_gather/blocksworld'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_blocksworld_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)
    
    save_input_dir = 'results_gather/BoxLift'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_boxlift_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)
    
    save_input_dir = 'results_gather/number_multiply'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_number_multiply_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)

    save_input_dir = 'results_gather/Letters'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_letters_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)
    
    save_input_dir = 'results_gather/game24'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_game24_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)

    save_input_dir = 'results_gather/path_plan'
    if not os.path.exists(save_input_dir):
        os.makedirs(save_input_dir)
    run_path_plan_baselines(save_input_dir, gather_save_input_dir, model_name, baseline_method_name, args_path)
    '''