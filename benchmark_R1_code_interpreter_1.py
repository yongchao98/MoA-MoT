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

from benchmark.run_R1_code_interpreter_data_syn import run_logic_game_baselines

if __name__ == '__main__':
    # gpt-4o, gpt-4o-mini, gpt-3.5-turbo for OpenAi API

    def log_run_info(log_file, run_info):
        with open(log_file, 'a') as f:
            f.write(run_info + "\n")

    infer_base_path = '/proj/long-multi/kqian/speech_2/LLaMA-Factory/examples/inference/'
    #infer_base_path = '/n/vlassak_lab/Lab/simulation_data/NLP_robotics/experiment/T5/large_model/llama3/LLaMA-Factory/examples/inference/'

    LLM_name_path_list = [['Qwen2.5-14B-1M', 'qwen2_14B.yaml'], ['Qwen2.5-7B-1M', 'qwen2_7B.yaml'], ['DeepSeek-R1-7B', 'deepseek_r1_7B.yaml'],
                          ['DeepSeek-R1-8B', 'deepseek_r1_8B.yaml'], ['Llama3-8B', 'llama3_8B.yaml'], ['gpt-4o', ''],
                          ['qwen2_7B_R1_CI_1', 'qwen2_7B_R1_CI_1.yaml'], ['llama3_8B_R1_CI_1', 'llama3_8B_R1_CI_1.yaml'],
                          ['qwen2_7B_R1_CI_2', 'qwen2_7B_R1_CI_2.yaml'], ['llama3_8B_R1_CI_2', 'llama3_8B_R1_CI_2.yaml']]

    for model_name, LLM_path in [['claude-3-5-sonnet-20241022', '']]:
        #model_name = 'gpt-4o'  # Qwen2.5-14B-1M, Qwen2.5-7B-1M, DeepSeek-R1-7B, DeepSeek-R1-8B, Llama3-8B
        # o3-mini-2025-01-31, gpt-4o, gpt-3.5-turbo, "claude-3-5-sonnet-20241022", o1, o1-preview, gpt-4o, qwen2_QwQ_32B-preview, 'DeepSeek-R1'

        args_path = os.path.join(infer_base_path, LLM_path)

        gather_save_input_dir = 'results_gather'
        baseline_method_name = '1_only_ques' # 1_only_ques, All_code_CoT, All_text,
        # R1_code_interpreter_test, R1_code_interpreter_data_syn_1, R1_code_interpreter_data_syn_hint, CodeSteer
        start_index = 20
        max_sample_num = 30
        max_round_num = 3
        ### Test baseline methods
        # filtered tasks: string_synthesis, letter_logic_diagram

        ### 33 tasks from SymBench
        env_name_list1 = ['eight_queens', 'pooling', 'reversi', 'light_puzzles', 'new_operator', 'mahjong_pattern', 'statistical_counting', 'synthesis_decomposition', '2048', 'matrix_transformation']

        env_name_list2 = ['pattern_recognition', 'constrained_linear_arrangement', 'string_synthesis', 'logic_puzzle', 'string_insertion', 'letter_logic_diagram', 'standard_sudoku', 'string_deletion_and_modification', 'string_splitting', 'cryptanalysis']

        env_name_list3 = ['permutations_and_combinations', 'logical_equation', 'combinatorial_calculation', 'game24', 'letters', 'number_multiply', 'gsm', 'math_geometry', 'math_counting_and_probability', 'BoxNet_v2']

        env_name_list4 = ['BoxLift', 'Blocksworld', 'Gridworld']

        ### 27 tasks from BigBench Hard
        env_name_big_bench_hard = ['big_bench_hard:' + task for task in ['date_understanding', 'web_of_lies', 'disambiguation_qa', 'formal_fallacies', 'geometric_shapes',
                      'logical_deduction_seven_objects', 'navigate', 'dyck_languages', 'boolean_expressions', 'causal_judgement',
                      'hyperbaton', 'logical_deduction_five_objects', 'logical_deduction_three_objects', 'movie_recommendation',
                      'multistep_arithmetic_two', 'object_counting', 'penguins_in_a_table', 'word_sorting', 'tracking_shuffled_objects_three_objects',
                      'tracking_shuffled_objects_seven_objects','tracking_shuffled_objects_five_objects', 'temporal_sequences',
                      'sports_understanding', 'snarks', 'salient_translation_error_detection', 'ruin_names', 'reasoning_about_colored_objects']]

        ### 100 Reasoning Gym tasks
        reasoning_gym_datasets = ['reasoning_gym_' + task for task in ['ab', 'acre', 'advanced_geometry', 'aiw', 'arc_1d', 'arc_agi']]

        base_path = 'results_gather'
        runtime_list = []
        for task_name in ['cryptanalysis']:
            start_time = time.time()
            run_logic_game_baselines(task_name, gather_save_input_dir, model_name, baseline_method_name, args_path, start_index, max_sample_num, max_round_num)
            end_time = time.time()
            runtime = end_time - start_time
            runtime_list.append(runtime/max_sample_num)

            output_path = base_path + f'/Cost_runtime_gather_{model_name}_{baseline_method_name}.txt'
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'a') as f:
                f.write(f"{runtime/max_sample_num}\n")
            print(f'Mean time cost: {np.mean(runtime_list)}')
            print(f'\nDataset saved to: {output_path}')