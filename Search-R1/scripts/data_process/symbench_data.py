"""
Preprocess the SymBench dataset to parquet format with task counts
"""

import re
import os
from r1_code_inter.Logic_Game_func import load_task_dataset
import argparse
import sys
import random
from datasets import Dataset
from collections import Counter
from typing import Union, List, Tuple, Dict, Optional

R1_code_interpreter_data_syn_prompt1 = r'''
The User asks a question, and you solve it. 
You first generate the reasoning and thinking process and then provide the User with the final answer.
During the thinking process, **you can generate python code** for efficient searching, optimization, and computing with the format of starting the python block with ```python. 
**A code query must involve only a single script that uses 'print' function for the output.**. 
Once the code script is complete, stop the generation. Then, the code interpreter platform will execute the code and return the execution output and error.
Once you feel you are ready for the final answer, directly return the answer with the format <<<answer content>>> at the end of your response.
Otherwise, you can continue your reasoning process and possibly generate more code query to solve the problem.

    '''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--local_dir', default='./data/nq_search')
    args = parser.parse_args()

    ### 33 tasks from SymBench
    env_name_list1 = ['pooling', 'reversi', 'light_puzzles', 'new_operator', 'mahjong_pattern',
                      'statistical_counting', 'synthesis_decomposition', '2048', 'matrix_transformation']
    env_name_list2 = ['pattern_recognition', 'constrained_linear_arrangement', 'string_synthesis', 'logic_puzzle',
                      'string_insertion', 'letter_logic_diagram', 'standard_sudoku', 'string_deletion_and_modification']
    env_name_list3 = ['string_splitting', 'permutations_and_combinations', 'logical_equation',
                      'combinatorial_calculation', 'cryptanalysis', 'BoxLift', 'Blocksworld', 'gsm', 'math_geometry']
    env_name_list4 = ['eight_queens', 'game24', 'letters', 'number_multiply', 'math_counting_and_probability', 'BoxNet_v2', 'Gridworld']

    ### 27 tasks from BigBench Hard
    env_name_big_bench_hard = ['big_bench_hard:' + task for task in
                               ['date_understanding', 'web_of_lies', 'disambiguation_qa', 'formal_fallacies',
                                'geometric_shapes',
                                'logical_deduction_seven_objects', 'navigate', 'dyck_languages', 'boolean_expressions',
                                'causal_judgement',
                                'hyperbaton', 'logical_deduction_five_objects', 'logical_deduction_three_objects',
                                'movie_recommendation',
                                'multistep_arithmetic_two', 'object_counting', 'penguins_in_a_table', 'word_sorting',
                                'tracking_shuffled_objects_three_objects',
                                'tracking_shuffled_objects_seven_objects', 'tracking_shuffled_objects_five_objects',
                                'temporal_sequences',
                                'sports_understanding', 'snarks', 'salient_translation_error_detection', 'ruin_names',
                                'reasoning_about_colored_objects']]

    ### 100 Reasoning Gym tasks
    reasoning_gym_available_datasets = [
        'ab', 'acre', 'advanced_geometry', 'aiw', 'arc_1d', 'arc_agi', 'base_conversion', 'basic_arithmetic', 'bf',
        'binary_alternation', 'binary_matrix', 'bitwise_arithmetic', 'caesar_cipher', 'calendar_arithmetic',
        'chain_sum', 'circuit_logic',
    ]
    reasoning_gym_datasets = ['reasoning_gym_' + task for task in reasoning_gym_available_datasets]

    # Counters to track task counts
    train_task_counter = Counter()
    test_task_counter = Counter()

    train_dataset_list = [];
    test_dataset_list = []
    idx_train = 0;
    idx_test = 0

    #for task_name in env_name_list1 + env_name_list2 + env_name_list4 + env_name_big_bench_hard + reasoning_gym_datasets:
    for task_name in env_name_list1 + env_name_list2 + env_name_list3 + env_name_list4 + env_name_big_bench_hard + reasoning_gym_datasets:
        solution_list, question_list, target_list, puzzles, solution_data_list, question_constrained_list, question_matrix_list, number_list, word_list, letter_list, save_input_dir = load_task_dataset(
            task_name, '')

        question_num_total = len(question_list)

        # Count for this specific task
        task_train_count = 0
        task_test_count = 0

        for i in range(min(5, question_num_total)):
            solution = '';
            number_list_item = '';
            target = '';
            word = '';
            letter = '';
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

            question = R1_code_interpreter_data_syn_prompt1 + 'question: ' + question + '\n'

            ground_truth = {"index": str(i), 'task_name': task_name}

            if task_name in env_name_list1 + env_name_list2 + env_name_list3 + env_name_big_bench_hard[
                                                                               :20] + reasoning_gym_datasets[:10]:
                idx = idx_train
                idx_train += 1
                split = 'train'
                task_train_count += 1
            elif task_name in env_name_list4 + env_name_big_bench_hard[20:] + reasoning_gym_datasets[10:]:
                idx = idx_test
                idx_test += 1
                split = 'test'
                task_test_count += 1

            data = {
                "data_source": task_name,
                "prompt": [{
                    "role": "user",
                    "content": question,
                }],
                "ability": "fact-reasoning",
                "reward_model": {
                    "style": "rule",
                    "ground_truth": ground_truth
                },
                "extra_info": {
                    'split': split,
                    'index': idx,
                }
            }

            if split == 'train':
                train_dataset_list.append(data)
            else:
                test_dataset_list.append(data)

        # Update counters with counts for this task
        if task_train_count > 0:
            train_task_counter[task_name] = task_train_count
        if task_test_count > 0:
            test_task_counter[task_name] = task_test_count

    # Print counts for each task
    print("\n===== TRAIN DATASET COUNTS =====")
    for task_name, count in sorted(train_task_counter.items()):
        print(f"{task_name}: {count} datapoints")
    print(f"Total train datapoints: {sum(train_task_counter.values())}")

    print("\n===== TEST DATASET COUNTS =====")
    for task_name, count in sorted(test_task_counter.items()):
        print(f"{task_name}: {count} datapoints")
    print(f"Total test datapoints: {sum(test_task_counter.values())}")

    # Save datasets to parquet files
    local_dir = args.local_dir
    os.makedirs(local_dir, exist_ok=True)

    random.shuffle(train_dataset_list)
    random.shuffle(test_dataset_list)

    train_ds = Dataset.from_list(train_dataset_list)  # Changed from train_list to train_dataset_list
    test_ds = Dataset.from_list(test_dataset_list)  # Changed from test_list to test_dataset_list
    train_ds.to_parquet(os.path.join(local_dir, 'train_symbench_3_simple.parquet'))
    test_ds.to_parquet(os.path.join(local_dir, 'test_symbench_3_simple.parquet'))

    '''
    # Save datasets as JSON Lines
    import json

    local_dir = args.local_dir
    os.makedirs(local_dir, exist_ok=True)

    train_path = os.path.join(local_dir, 'train_symbench.jsonl')
    test_path  = os.path.join(local_dir, 'test_symbench.jsonl')

    with open(train_path, 'w', encoding='utf-8') as f:
        for rec in train_dataset_list:
            json.dump(rec, f, ensure_ascii=False)
            f.write('\n')

    with open(test_path, 'w', encoding='utf-8') as f:
        for rec in test_dataset_list:
            json.dump(rec, f, ensure_ascii=False)
            f.write('\n')

    print(f"\nDatasets saved to {local_dir}")
    '''
