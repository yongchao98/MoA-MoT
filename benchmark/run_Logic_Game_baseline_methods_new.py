import json
import re
import pandas as pd
import os
import subprocess
import sys
from openai import OpenAI
from generation_models import message_construct_func, message_construct_llama_func, GPT_response, count_total_tokens, \
    extract_code, extract_and_check, LLM_answer_code_checker, paraphrase_with_GPT4, log_run_info, load_conversation_data
import random
import math
import json
from typing import List, Tuple, Dict
import time
import numpy as np
import ast
from prompt import *
from argparse import ArgumentParser
from symbolic_code_check import analyze_computational_approach, analyze_code_and_explain
import random
import string
from LLaMA_Factory.src.llamafactory.chat.chat_model import run_response

### To do, add tasks to import related functions
from Logic_Game_func import load_task_dataset, verify_solution_func_gather, multi_round_answer_sampling

def save_file_func_baselines(save_code_dir, response_list, user_prompt_list, question, system_message):
    data = {
        'question': question,
        'response_list': response_list,
        'user_prompt_list': user_prompt_list,
        'system_message': system_message
    }

    output_file = os.path.join(save_code_dir, 'conversation_data.json')

    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"Data successfully saved to {output_file}")
    except Exception as e:
        print(f"Error saving data: {str(e)}")

def run_logic_game_baselines(task_name, gather_save_input_dir, model_name, baseline_method_name, args_path, max_sample_num):
    print('\n' + '*'*30)
    total_sample_num = 0
    total_correct_num = 0

    # Load dataset
    solution_list, question_list, target_list, puzzles, solution_data_list, question_constrained_list, question_matrix_list, number_list, word_list, letter_list, save_input_dir = load_task_dataset(
        task_name, model_name)

    print(f'length of number_list: {len(number_list)}')
    base_save_code_dir = save_input_dir + f'/result_{task_name}_{baseline_method_name}_{model_name}'
    if not os.path.exists(base_save_code_dir):
        os.makedirs(base_save_code_dir)

    ### Remain unchanged
    for i in range(0, min(max(len(solution_list), len(question_list), len(number_list)), max_sample_num)):
        solution = ''; question = ''; number_list_item = ''; target = ''; word = ''; letter = ''
        question = question_list[i]
        if len(solution_list) > 0:
            solution = solution_list[i]
        if len(number_list) > 0:
            number_list_item = number_list[i]
        if len(target_list) > 0:
            target = target_list[i]
        if len(word_list) > 0:
            word = word_list[i]
        if len(letter_list) > 0:
            letter = letter_list[i]

        total_sample_num += 1

        save_code_dir = os.path.join(base_save_code_dir, f"Test_sample_{i}/")
        if not os.path.exists(save_code_dir):
            os.makedirs(save_code_dir)

        print('-------###-------###-------###-------')
        print(f"\nTest_sample_{i}, total: {min(len(solution_list), max_sample_num)}/")

        code_interpreter = False
        system_message = ''
        if baseline_method_name == '1_only_ques':
            user_prompt_list = [question]
        elif baseline_method_name == 'code_interpreter':
            code_interpreter = True
            user_prompt_list = [question]
        elif baseline_method_name == 'AutoGen':
            user_prompt_list = [AutoGen_prompt + question]
        elif baseline_method_name == 'All_code_CoT':
            user_prompt_list = [with_COT_code_output_prompt + question]
        elif baseline_method_name == 'All_text':
            user_prompt_list = [text_output_prompt + question]

        if model_name in ['o3-mini-2025-01-31', 'o1', "o1-preview", 'o1-mini', 'gpt-4o', 'gpt-4o-mini', 'gpt-3.5-turbo', "claude-3-5-sonnet-20241022", "claude-3-sonnet-20240229",
                             "claude-3-opus-20240229", "claude-3-haiku-20240307", 'open-mixtral-8x7b', "mistral-large-latest", 'DeepSeek-R1']:
            response = GPT_response(system_message, user_prompt_list[0], model_name=model_name, code_interpreter=code_interpreter,
                                         user_prompt_list=user_prompt_list, response_total_list=[], logprobs=False)
        else:
            messages = message_construct_llama_func(user_prompt_list, [])
            response = run_response(messages, args_path)
        response_list = []
        response_list.append(response)

        save_file_func_baselines(save_code_dir, response_list, user_prompt_list, question, system_message)

        response = response_list[-1]
        if count_total_tokens([response], []) > 10000:
            response = response[:10000]

        True_false_result_1, True_false_result_2 = verify_solution_func_gather(i, task_name, response,
                                                                               save_code_dir, question, solution,
                                                                               target, puzzles, solution_data_list,
                                                                               solution_list, question_constrained_list,
                                                                               question_matrix_list, number_list_item, word, letter)

        ### Remain unchanged
        if True_false_result_1 == False and True_false_result_2 == False:
            print('False')
            with open(save_code_dir + f"/success_failure.txt", "w") as f:
                f.write('False')
        else:
            print('True')
            with open(save_code_dir + f"/success_failure.txt", "w") as f:
                f.write('True')
            total_correct_num += 1

        print(f'\ntotal_sample_num: {total_sample_num}')
        print(f'total_correct_num: {total_correct_num}\n')

        with open(base_save_code_dir + f"/acc_result_log_{model_name}_{baseline_method_name}.txt", "w") as f:
            f.write(f"correct/all:{total_correct_num}/{total_sample_num}\n")

    run_info = f"CodeSteer, {task_name}, {baseline_method_name}, {model_name}\n"
    run_info_result = f'correct/all:{total_correct_num}/{total_sample_num}\n'
    log_file_result = os.path.join(gather_save_input_dir, f"acc_result_log_{model_name}.txt")
    log_run_info(log_file_result, run_info + run_info_result)