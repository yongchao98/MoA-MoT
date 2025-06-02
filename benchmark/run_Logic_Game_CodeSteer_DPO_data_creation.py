import json
import re
import pandas as pd
import os
import subprocess
import sys
from openai import OpenAI
from generation_models import message_construct_func, message_construct_llama_func, GPT_response, count_total_tokens, extract_code, \
    extract_and_check, LLM_answer_code_checker, save_file_func, paraphrase_with_GPT4, log_run_info, load_conversation_data, save_dataset_item
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
import copy
#from LLaMA_Factory.src.llamafactory.chat.chat_model import run_response

### To do, add tasks to import related functions
from Logic_Game_func import load_task_dataset, verify_solution_func_gather, multi_round_answer_sampling

def run_logic_game_DPO_data(task_name, gather_save_input_dir, model_name, max_tree_depth, args_path, args_path_DPO_guidance_prompt_gen, CodeSteer_LLM, CodeSteer_LLM_2, max_sample_num, DPO_sampling = False, original_answer_base_save_code_dir = None):
    print('\n' + '*'*30)
    total_sample_num = 0
    total_correct_num = 0

    # Load dataset
    solution_list, question_list, target_list, puzzles, solution_data_list, question_constrained_list, question_matrix_list = load_task_dataset(
        task_name, model_name)

    base_save_code_dir = save_input_dir + f'/result_logical_game_{CodeSteer_LLM_1}_{CodeSteer_LLM_2}_{model_name}_MTD_{max_tree_depth}_CodeSteer_DPO_1'

    ### Remain unchanged
    for i in range(0, min(len(solution_list), max_sample_num)):
        question = question_list[i]
        solution = solution_list[i]
        target = target_list[i]
        total_sample_num += 1
        if puzzles == []:
            puzzle = None
        else:
            puzzle = puzzles[i]

        if original_answer_base_save_code_dir:
            print(f'\nUse original answer base save code dir: {original_answer_base_save_code_dir}\n')
            base_save_code_dir = original_answer_base_save_code_dir
            print(base_save_code_dir)
            print(f"Test_sample_{i}/")
            save_code_dir = os.path.join(base_save_code_dir, f"Test_sample_{i}/")
            if not os.path.exists(save_code_dir):
                print(f'The original save_code_dir does not exist: {save_code_dir}')
                continue

        else:
            if not os.path.exists(base_save_code_dir):
                os.makedirs(base_save_code_dir)
            print(f'\nUse new answer base save code dir: {base_save_code_dir}\n')
            print(base_save_code_dir)
            print(f"Test_sample_{i}/")
            save_code_dir = os.path.join(base_save_code_dir, f"Test_sample_{i}/")
            if not os.path.exists(save_code_dir):
                os.makedirs(save_code_dir)

            print('-------###-------###-------###-------')
            print(f"\nTest_sample_{i}, total: {len(solution_list)}/")
            response_list = multi_round_answer_sampling(save_code_dir, question, [], [], [], [], [], model_name, CodeSteer_LLM, args_path, max_tree_depth, 0)
            True_false_result_1, True_false_result_2 = verify_solution_func_gather(i, task_name, response_list[-1], save_code_dir, question, solution, target, puzzles, solution_data_list, solution_list, question_constrained_list, question_matrix_list)

            if True_false_result_1 or True_false_result_2:
                total_correct_num += 1
            print(f'\ntotal_sample_num: {total_sample_num}')
            print(f'total_correct_num: {total_correct_num}\n')
            with open(base_save_code_dir + f"/acc_result_log_{model_name}.txt", "w") as f:
                f.write(f"correct/all:{total_correct_num}/{total_sample_num}\n")

        if DPO_sampling:
            ## DPO sampling
            conv_path = os.path.join(save_code_dir, 'conversation_data.json')
            status_path = os.path.join(save_code_dir, 'success_failure.txt')
            response_list, user_prompt_list, question, CodeSteer_input_prompt_list, \
            CodeSteer_input_prompt_training_list, CodeSteer_output_prompt_guidance_list = \
                load_conversation_data(conv_path)
            with open(status_path, 'r') as f:
                status = f.read().strip()

            if status == 'False':
                score_original = -len(CodeSteer_output_prompt_guidance_list)
            else:
                score_original = 15 - len(CodeSteer_output_prompt_guidance_list)

            Get_sample_pair = False
            for current_step_num in range(0, len(response_list)):
                if Get_sample_pair:
                    break
                for iteration_index in range(min(5 - current_step_num, 2)):
                    print(f'\nTotal length of response_list: {len(response_list)}\n')
                    print(f'DPO sample: {i}_{current_step_num}_{iteration_index}')
                    save_code_dir_DPO_sample = os.path.join(base_save_code_dir, f"Test_sample_{CodeSteer_LLM_2}_{i}_{current_step_num}_{iteration_index}/")
                    if not os.path.exists(save_code_dir_DPO_sample):
                        os.makedirs(save_code_dir_DPO_sample)
                    response_list_sampled = multi_round_answer_sampling(save_code_dir_DPO_sample, question, response_list[:current_step_num],
                                                                CodeSteer_output_prompt_guidance_list[:current_step_num],
                                                                CodeSteer_input_prompt_list[:current_step_num],
                                                                CodeSteer_input_prompt_training_list[:current_step_num], user_prompt_list[:current_step_num], model_name,
                                                                CodeSteer_LLM_2, args_path_DPO_guidance_prompt_gen, max_tree_depth,
                                                                current_step_num)
                    True_false_result_1, True_false_result_2 = verify_solution_func_gather(i, task_name, response_list_sampled[-1], save_code_dir_DPO_sample, question, solution, target, puzzles, solution_data_list, solution_list, question_constrained_list, question_matrix_list)

                    conv_path = os.path.join(save_code_dir_DPO_sample, 'conversation_data.json')
                    status_path = os.path.join(save_code_dir_DPO_sample, 'success_failure.txt')
                    response_list_new, user_prompt_list_new, question, CodeSteer_input_prompt_list_new, \
                    CodeSteer_input_prompt_training_list_new, CodeSteer_output_prompt_guidance_list_new = \
                        load_conversation_data(conv_path)
                    with open(status_path, 'r') as f:
                        status_new = f.read().strip()
                    if status_new == 'False':
                        score_new = -len(CodeSteer_output_prompt_guidance_list)
                    else:
                        score_new = 15 - len(CodeSteer_output_prompt_guidance_list)

                    print(f'\n{i}_{current_step_num}_{iteration_index}: score_original:{score_original} -> score_new:{score_new}\n')

                    if score_new != score_original:
                        max_token_len = 7000
                        if count_total_tokens(CodeSteer_input_prompt_training_list[:current_step_num],
                                              CodeSteer_output_prompt_guidance_list[:current_step_num]) > max_token_len:
                            continue
                        # Verify we have matching data
                        if len(CodeSteer_input_prompt_training_list) == len(CodeSteer_output_prompt_guidance_list):
                            # Create history list from previous interactions
                            conversations = []
                            for index in range(len(CodeSteer_input_prompt_training_list[:current_step_num])):
                                conversations.append({"from": "human", "value": CodeSteer_input_prompt_training_list[index]})
                                conversations.append({"from": "gpt", "value": CodeSteer_output_prompt_guidance_list[index]})
                            conversations.append(
                                {"from": "human", "value": CodeSteer_input_prompt_training_list[current_step_num]})

                            if score_new > score_original and score_new > 0:
                                dataset_item = {
                                    "conversations": conversations,
                                    "chosen": {"from": "gpt", "value": CodeSteer_output_prompt_guidance_list_new[current_step_num]},
                                    "rejected": {"from": "gpt", "value": CodeSteer_output_prompt_guidance_list[current_step_num]},
                                    "score_chosen": score_new,
                                    "score_rejected": score_original
                                }
                                save_dataset_item(dataset_item, gather_save_input_dir + f"/CodeSteer_DPO_dataset_{task_name}_1.json")
                                Get_sample_pair = True
                                break
                            elif score_new < score_original and score_original > 0:
                                dataset_item = {
                                    "conversations": conversations,
                                    "chosen": {"from": "gpt", "value": CodeSteer_output_prompt_guidance_list[current_step_num]},
                                    "rejected": {"from": "gpt", "value": CodeSteer_output_prompt_guidance_list_new[current_step_num]},
                                    "score_chosen": score_original,
                                    "score_rejected": score_new
                                }
                                save_dataset_item(dataset_item, gather_save_input_dir + f"/CodeSteer_DPO_dataset_{task_name}_1.json")
                                Get_sample_pair = True
                                break

    run_info = f"CodeSteer, {task_name}, {CodeSteer_LLM}, {CodeSteer_LLM_2}, {model_name}, MTD_{max_tree_depth}_CodeSteer_1\n"
    run_info_result = f'correct/all:{total_correct_num}/{total_sample_num}\n'
    log_file_result = os.path.join(gather_save_input_dir, f"acc_result_log_{model_name}.txt")
    log_run_info(log_file_result, run_info + run_info_result)