import os
import json
from typing import Union, List, Tuple, Dict, Optional
import pandas as pd
import subprocess
import sys
from openai import OpenAI
from generation_models import message_construct_func, message_construct_llama_func, GPT_response, count_total_tokens, extract_code, extract_and_check, LLM_answer_code_checker, save_file_func, paraphrase_with_GPT4, log_run_info, load_conversation_data
import random
from typing import List, Tuple, Dict
import time
import numpy as np
import ast
import string
import copy
import re
import math
from dataclasses import dataclass
from prompt import *
from symbolic_code_check import analyze_computational_approach, analyze_code_and_explain
#from LLaMA_Factory.src.llamafactory.chat.chat_model import run_response

##### Logic Game #####
### All the related functions for each task are listed in order below.

##### Gathered functions for all tasks #####
def multi_round_answer_sampling(save_code_dir, question, response_list_input, CodeSteer_output_prompt_guidance_list_input, CodeSteer_input_prompt_list_input, CodeSteer_input_prompt_training_list_input,
                                user_prompt_list_input, model_name, CodeSteer_LLM, args_path, max_tree_depth, current_tree_depth):
    print(f'\nLength of user_prompt_list: {len(user_prompt_list_input)}')
    print(f'Length of response_list: {len(response_list_input)}\n')

    if current_tree_depth == 0:
        response_list_input = [];
        CodeSteer_output_prompt_guidance_list_input = [];
        CodeSteer_input_prompt_list_input = [code_text_choice_prompt + question];
        CodeSteer_input_prompt_training_list_input = [code_text_choice_prompt + question]

        ############ Starting first guidance ############
        if CodeSteer_LLM.startswith("llama"):
            print(f'\nCodeSteer_LLM: {CodeSteer_LLM}, open model!\n')
            messages = message_construct_llama_func([code_text_choice_prompt + question], [])
            starting_prompt_choice = run_response(messages, args_path)
            print(f'Starting prompt choice: {starting_prompt_choice}')
            user_prompt_list_input = [starting_prompt_choice + question]
            CodeSteer_output_prompt_guidance_list_input.append(starting_prompt_choice)
        elif CodeSteer_LLM in ['gpt-4o', 'gpt-4o-mini', 'gpt-3.5-turbo', "claude-3-sonnet-20240229",
                               "claude-3-opus-20240229", "claude-3-haiku-20240307"]:
            print('\nAPI based close models!\n')
            starting_prompt_choice = GPT_response("", code_text_choice_prompt + question, model_name=CodeSteer_LLM,
                                                  code_interpreter=False,
                                                  user_prompt_list=[code_text_choice_prompt + question],
                                                  response_total_list=[], logprobs=False, temperature=1.0)
            print(f'Starting prompt choice: {starting_prompt_choice}')
            if '<<<Code>>>' in starting_prompt_choice:
                # paraphrased_prompt_guidance = paraphrase_with_GPT4(with_COT_code_output_prompt)
                paraphrased_prompt_guidance = with_COT_code_output_prompt
                user_prompt_list_input = [paraphrased_prompt_guidance + question]
            elif '<<<Text>>>' in starting_prompt_choice:
                # paraphrased_prompt_guidance = paraphrase_with_GPT4(text_output_prompt)
                paraphrased_prompt_guidance = text_output_prompt
                user_prompt_list_input = [paraphrased_prompt_guidance + question]
            else:
                print('Error: No text/code choice made')
                raise ValueError('Error: No text/code choice made')
            CodeSteer_output_prompt_guidance_list_input.append(paraphrased_prompt_guidance)
        else:
            raise ValueError('Error: No LLM model specified')

        response = GPT_response('', user_prompt_list_input[0], model_name=model_name, code_interpreter=False,
                                user_prompt_list=user_prompt_list_input, response_total_list=response_list_input, logprobs=False)
        response_list_input.append(response)
        # print(f'\nResponse_0: {response}\n')

    ############ Further rounds of guidance ############
    for tree_depth in range(max_tree_depth):
        response = response_list_input[-1]
        code_block_list = extract_code(response)
        if len(code_block_list) > 0:
            code_complexity_summary, code_complexity_score = analyze_code_and_explain(code_block_list[0])
            if code_complexity_score <= 2:
                code_complexity_summary += '\nThe generated code may not be complex enough to carry out symbolic computing for solving the task.'
            with open(save_code_dir + f"/code_1_{tree_depth}.py", "w") as f:
                f.write(code_block_list[0])

            try:
                result = subprocess.run(
                    ["python3", save_code_dir + f"/code_1_{tree_depth}.py"],
                    capture_output=True, text=True, timeout=10
                )
                output = result.stdout
                errors = result.stderr
            except subprocess.TimeoutExpired as e:
                output = e.stdout if e.stdout else ""
                errors = e.stderr if e.stderr else ""
                errors += f"\nTimeoutExpired: Command '{e.cmd}' timed out after {e.timeout} seconds"

            if isinstance(output, str):
                if count_total_tokens([output + errors], []) > 2000:
                    response = response + f'\nThe execution result from the generated code is too long to be displayed.'
                else:
                    response = response + f'\nThe execution result from the generated code is:\noutput: {output}, errors: {errors}'
            else:
                response = response + f'\nThe execution result from the generated code is:\nerrors: {errors}'

        check_code_saving_path = save_code_dir + f"/check_code_1_{tree_depth}.py"
        check_result = LLM_answer_code_checker(question, response, check_code_saving_path)

        CodeSteer_input_prompt_head = f'''{decision_prompt} {question}\n'''
        if len(code_block_list) > 0:
            print('\n############True#############\n')
            CodeSteer_input_prompt = f'''The response from TaskLLM is: {response}\n\nThe feedback from the checking agent is:\n{check_result}\n\nThe summary of generated code complexity is: {code_complexity_summary}\n\n''' + \
                                     f'''The final returned guidance prompt should be of the format <<<guidance prompt content>>>.'''

        else:
            CodeSteer_input_prompt = f'''The response from TaskLLM is: {response}\n\nThe feedback from the checking agent is:\n{check_result}\n\n''' + \
                                     f'''The final returned guidance prompt should be of the format <<<guidance prompt content>>>.'''

        CodeSteer_input_prompt_total = CodeSteer_input_prompt_head + CodeSteer_input_prompt
        CodeSteer_input_prompt_list_input.append(CodeSteer_input_prompt_total)
        CodeSteer_input_prompt_training_list_input.append(CodeSteer_input_prompt)

        if CodeSteer_LLM.startswith("llama"):
            print(f'\nCodeSteer_LLM: {CodeSteer_LLM}, open model!\n')
            messages = message_construct_llama_func(CodeSteer_input_prompt_training_list_input,
                                                    CodeSteer_output_prompt_guidance_list_input)
            guidance_prompt = run_response(messages, args_path)
        elif CodeSteer_LLM in ['gpt-4o', 'gpt-4o-mini', 'gpt-3.5-turbo', "claude-3-sonnet-20240229",
                               "claude-3-opus-20240229", "claude-3-haiku-20240307"]:
            print('\nAPI based close models!\n')
            response_text = GPT_response("", '', model_name=CodeSteer_LLM, code_interpreter=False,
                                         user_prompt_list=CodeSteer_input_prompt_list_input,
                                         response_total_list=CodeSteer_output_prompt_guidance_list_input,
                                         logprobs=False, temperature=1.0)

            matches = re.findall(r'<<<(.*?)>>>', response_text, re.DOTALL)
            guidance_prompt = matches[-1] if matches else response_text
        else:
            raise ValueError('Error: No LLM model specified')
        print(f'\nGuidance prompt_{tree_depth + 1}: {guidance_prompt}\n')

        CodeSteer_output_prompt_guidance_list_input.append(guidance_prompt)
        if '<<<Code>>>' in guidance_prompt:
            guidance_prompt = with_COT_code_output_prompt
        elif '<<<Text>>>' in guidance_prompt:
            guidance_prompt = text_output_prompt
        elif '<<<Return Answer>>>' in guidance_prompt or 'Return Answer' in guidance_prompt or '<<<Terminate>>>' in guidance_prompt or 'Terminate' in guidance_prompt:
            break
        user_prompt_list_input.append(guidance_prompt)

        response = GPT_response('', user_prompt_list_input[0], model_name=model_name, code_interpreter=False,
                                user_prompt_list=user_prompt_list_input, response_total_list=response_list_input, logprobs=False)

        response_list_input.append(response)
        # print(f'\nResponse_{tree_depth}: {response}\n')
    save_file_func(save_code_dir, response_list_input, user_prompt_list_input, question, CodeSteer_input_prompt_list_input,
                   CodeSteer_input_prompt_training_list_input, CodeSteer_output_prompt_guidance_list_input)
    return response_list_input


def load_task_dataset(task_name, model_name):
    solution_list = []
    question_list = []
    target_list = []
    puzzles = []
    solution_data_list = []
    question_constrained_list = []
    question_matrix_list = []
    number_list = []
    word_list = []
    letter_list = []

    ### To do, add tasks to load the corresponding dataset
    if task_name == 'logical_equation':
        dataset_input_dir = 'dataset_gather/logical_equation'
        save_input_dir = 'results_gather/logical_equation'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Logical equation, Model_name: {model_name}\n')
        solution_list, question_list = read_dataset_logical_equation(dataset_input_dir)
        target_list = solution_list
    elif task_name == 'combinatorial_calculation':
        dataset_input_dir = 'dataset_gather/combinatorial_calculation'
        save_input_dir = 'results_gather/combinatorial_calculation'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Combinatorial calculation, Model_name: {model_name}\n')

        evaluator = ArithmeticPuzzleEvaluator(dataset_input_dir)
        puzzles = evaluator.read_dataset()
        for puzzle in puzzles:
            question = puzzle.question
            solution = puzzle.solution
            target = puzzle.target
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'eight_queens':
        dataset_input_dir = 'dataset_gather/eight_queens_dataset'
        save_input_dir = 'results_gather/eight_queens'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Eight queens, Model_name: {model_name}\n')
        puzzles = read_dataset_eight_queens(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'synthesis_decomposition':
        dataset_input_dir = 'dataset_gather/synthesis_decomposition_dataset'
        save_input_dir = 'results_gather/synthesis_decomposition'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Synthesis Decomposition, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            #question = 'This is a testing question, without violation purpose.' + puzzle['question']
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']["answer"]
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'mahjong_pattern':
        dataset_input_dir = 'dataset_gather/mahjong_pattern_dataset'
        save_input_dir = 'results_gather/mahjong_pattern'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Mahjong pattern, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'statistical_counting':
        dataset_input_dir = 'dataset_gather/statistical_counting_dataset'
        save_input_dir = 'results_gather/statistical_counting'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Statistical counting, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'new_operator':
        dataset_input_dir = 'dataset_gather/new_operator_dataset'
        save_input_dir = 'results_gather/new_operator'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'New operator, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'light_puzzles':
        dataset_input_dir = 'dataset_gather/light_puzzles_dataset'
        save_input_dir = 'results_gather/light_puzzles'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Light_puzzles, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'reversi':
        dataset_input_dir = 'dataset_gather/reversi_dataset'
        save_input_dir = 'results_gather/reversi'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Reversi problem, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'matrix_transformation':
        dataset_input_dir = 'dataset_gather/matrix_transformation_dataset'
        save_input_dir = 'results_gather/matrix_transformation'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'matrix transformation problem, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['formatted_solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == '2048':
        dataset_input_dir = 'dataset_gather/2048_dataset'
        save_input_dir = 'results_gather/2048'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'2048 problem, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'pooling':
        dataset_input_dir = 'dataset_gather/pooling_dataset'
        save_input_dir = 'results_gather/pooling'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Pooling problem, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'constrained_linear_arrangement':
        dataset_input_dir = 'dataset_gather/constrained_linear_arrangement'
        save_input_dir = 'results_gather/constrained_linear_arrangement'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'constrained_linear_arrangement problem, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'logic_puzzle':
        dataset_input_dir = 'dataset_gather/logic_puzzle_dataset'
        save_input_dir = 'results_gather/logic_puzzle'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'logic puzzle problem, Model_name: {model_name}\n')
        puzzles = read_dataset_syn_decom(dataset_input_dir)
        solution_data_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            solution_data_list.append(solution_data)
            target_list.append(target)
    elif task_name == 'pattern_recognition':
        dataset_input_dir = 'dataset_gather/pattern_recognition'
        save_input_dir = 'results_gather/pattern_recognition'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Pattern recognition, Model_name: {model_name}\n')
        puzzles = read_dataset_pattern_recognition(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            target_row, target_col = solution_data['row'], solution_data['column']
            solution_list.append([target_row, target_col])
            question_list.append(question)
            target_list = solution_list
    elif task_name == 'string_insertion':
        dataset_input_dir = 'dataset_gather/string_insertion'
        save_input_dir = 'results_gather/string_insertion'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'String Insertion, Model_name: {model_name}\n')
        puzzles = read_dataset_string_insertion(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['modified_string']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'letter_logic_diagram':
        dataset_input_dir = 'dataset_gather/letter_logic_diagram'
        save_input_dir = 'results_gather/letter_logic_diagram'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Letter_logic_diagram, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_letter_logic_diagram(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'string_synthesis':
        dataset_input_dir = 'dataset_gather/string_synthesis'
        save_input_dir = 'results_gather/string_synthesis'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'String Synthesis, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_string_synthesis(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['final_counts']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'standard_sudoku':
        dataset_input_dir = 'dataset_gather/standard_sudoku'
        save_input_dir = 'results_gather/standard_sudoku'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Standard Sudoku, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_standard_sudoku(dataset_input_dir)
        question_matrix_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['solution']
            question_matrix = solution_data['puzzle']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            question_matrix_list.append(question_matrix)
            target_list.append(target)
    elif task_name == 'permutations_and_combinations':
        dataset_input_dir = 'dataset_gather/permutations_and_combinations'
        save_input_dir = 'results_gather/permutations_and_combinations'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'permutations_and_combinations, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_permutations_and_combinations(dataset_input_dir)
        question_constrained_list = []
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['correct_solution']
            question_constrained_linear = solution_data['constraints_text']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            question_constrained_list.append(question_constrained_linear)
            target_list.append(target)
    elif task_name == 'string_deletion_and_modification':
        dataset_input_dir = 'dataset_gather/string_deletion_and_modification'
        save_input_dir = 'results_gather/string_deletion_and_modification'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'String Deletion And Modification, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_string_deletion_and_modification(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['final_string']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'minesweeper':
        dataset_input_dir = 'dataset_gather/minesweeper'
        save_input_dir = 'results_gather/minesweeper'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Minesweeper, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_minesweeper(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['mines']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'cryptanalysis':
        dataset_input_dir = 'dataset_gather/cryptanalysis'
        save_input_dir = 'results_gather/cryptanalysis'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'Cryptanalysis, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_cryptanalysis(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['answer']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'string_splitting':
        dataset_input_dir = 'dataset_gather/string_splitting'
        save_input_dir = 'results_gather/string_splitting'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        print(f'String Splitting, Model_name: {model_name}, CodeSteer\n')
        puzzles = read_dataset_string_splitting(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzle['question']
            solution_data = puzzle['solution_data']
            solution = solution_data['expected_answer']
            target = solution
            solution_list.append(solution)
            question_list.append(question)
            target_list.append(target)
    elif task_name == 'game24':
        dataset_input_dir = 'dataset_gather/game24_dataset/24'
        dataset_csv = dataset_input_dir + f'/24.csv'
        df = pd.read_csv(dataset_csv)
        question_prompt = f'Use numbers and basic arithmetic operations (+ - * /) to obtain 24. Each number should be used only once but each number has to be used in the equation. ' \
                          f'Input: 9 10 11 13, Answer: ((10-9)*(11+13)) = 24 Input: 4 10 10 11, Answer: ((4*11)-(10+10)) = 24 Input: 5 6 13 13, Answer: ((5-(13/13))*6)' \
                          f'Input: 2 6 6 7, Answer: ((6+(6*7))/2) Input: 2 6 10 18, Answer: (2-(6-(10+18)))'
        save_input_dir = 'results_gather/game24'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        for index_game24 in range(0, min(len(df), 1500), 5):
            number_list_item = list(map(int, df['Puzzles'][index_game24].split()))
            question = f'{question_prompt}'
            question += f'Input: '
            for number in number_list_item:
                question += f'{number} '
            question += f'Answer:\nOutput final answer with the format <<<answer>>>.'
            question_list.append(question)
            number_list.append(number_list_item)
    elif task_name == 'letters':
        dataset_input_dir = 'dataset_gather'
        save_input_dir = 'results_gather/letters'
        for min_length, max_length in [(10, 15), (15, 20), (20, 25)]:
            base_dir = dataset_input_dir + f'/Letters_dataset_min_length_{min_length}_max_length_{max_length}/'
            for i, letter in enumerate(string.ascii_lowercase[::6]):
                for letter_freq in range(1, 6):
                    for index in range(1):
                        saving_dir = base_dir + f"{letter}_{letter_freq}_{index}/"
                        word = read_words_from_file_letters(saving_dir + 'test_words.json')
                        question = create_prompt_letters(word, letter)
                        question_list.append(question)
                        word_list.append(word)
                        letter_list.append(letter)
    elif task_name == 'number_multiply':
        dataset_input_dir = 'dataset_gather/number_multiply'
        save_input_dir = 'results_gather/number_multiply'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        puzzles = read_dataset_number_multipy(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzles['question']
            question_list.append(question)
            solution = puzzle['solution_data']
            solution_list.append(solution)
    elif task_name == 'big_bench_hard':
        dataset_input_dir = 'dataset_gather/BIG-Bench-Hard/bbh'
        save_input_dir = 'results_gather/big_bench_hard'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
        puzzles = read_dataset_big_bench_hard(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzles['question']
            question_list.append(question)
            solution = puzzle['solution_data']
            solution_list.append(solution)
    elif task_name == 'gsm':
        dataset_input_dir = 'dataset_gather/gsmhardv2.jsonl'
        save_input_dir = 'results_gather/gsm'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
            puzzles = read_dataset_gsm(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzles['question']
            question_list.append(question)
            solution = puzzle['solution_data']
            solution_list.append(solution)
    elif task_name == 'math_geometry':
        dataset_input_dir = 'dataset_gather/MATH/train/geometry'
        save_input_dir = 'results_gather/math_geometry'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
            puzzles = read_dataset_math_geometry(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzles['question']
            question_list.append(question)
            solution = puzzle['solution_data']
            solution_list.append(solution)
    elif task_name == 'math_counting_and_probability':
        dataset_input_dir = 'dataset_gather/MATH/train/counting_and_probability'
        save_input_dir = 'results_gather/math_counting_and_probability'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
            puzzles = read_dataset_math_counting_and_probability(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzles['question']
            question_list.append(question)
            solution = puzzle['solution_data']
            solution_list.append(solution)
    elif task_name == 'BoxNet1':
        dataset_input_dir = 'dataset_gather/BoxNet1_dataset'
        save_input_dir = 'results_gather/BoxNet1'
        if not os.path.exists(save_input_dir):
            os.makedirs(save_input_dir)
            puzzles = read_dataset_boxnet1(dataset_input_dir)
        for puzzle in puzzles:
            question = puzzles['question']
            question_list.append(question)
            solution = puzzle['pg_dict']
            solution_list.append(solution)

    return solution_list, question_list, target_list, puzzles, solution_data_list, question_constrained_list, question_matrix_list, number_list, word_list, letter_list, save_input_dir

def verify_solution_func_gather(i, task_name, response, save_code_dir, question, solution, target, puzzles, solution_data_list, solution_list, question_constrained_list, question_matrix_list, number_list_item, word, letter):
    # Verify solution based on task type
    original_response = response

    code_block_list = extract_code(response)
    for index, code_string in enumerate(code_block_list):
        with open(save_code_dir + f"/code_1_{index}.py", "w") as f:
            f.write(code_string)

    # Test the generated code
    if not os.path.exists(save_code_dir + f"/code_1_0.py"):
        pass
    else:
        try:
            result = subprocess.run(
                ["python3", "-c", f"exec(open('{save_code_dir}/code_1_0.py').read()); print(result)"],
                capture_output=True, text=True, timeout=15)
            if result.stdout == '':
                result = subprocess.run(
                    ["python3", "-c", f"exec(open('{save_code_dir}/code_1_0.py').read()); print(Answer)"],
                    capture_output=True, text=True, timeout=15)
            if result.stdout == '':
                result = subprocess.run(
                    ["python3", "-c", f"exec(open('{save_code_dir}/code_1_0.py').read()); print(answer)"],
                    capture_output=True, text=True, timeout=15)

            # if '<<<' in result.stdout and '>>>' in result.stdout:
            response = result.stdout
            errors = result.stderr
        except Exception as e:
            pass

    ### To do, add tasks to extract the solution from the generated response and verify the correctness, here we check both response and original_response
    ### Take care the names of variables to be the same as other tasks since we need to save several variables in the following process.
    if task_name == 'logical_equation':
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_logical_equation(response)
        extracted_text_1, _ = extract_and_check(output_1)
        extracted_text_1 = extracted_text_1.strip()
        try:
            solution_1 = eval(extracted_text_1)
        except:
            solution_1 = []

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_logical_equation(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        extracted_text_2 = extracted_text_2.strip()
        try:
            solution_2 = eval(extracted_text_2)
        except:
            solution_2 = []

        constraints = [line.split('. ')[1] for line in question.split('\n')
                       if line.strip() and line[0].isdigit()]

        True_false_result_1 = verify_solution_logical_equation(solution, solution_1, constraints)
        True_false_result_2 = verify_solution_logical_equation(solution, solution_2, constraints)
    elif task_name == 'combinatorial_calculation':
        dataset_input_dir = 'dataset_gather/combinatorial_calculation'
        evaluator = ArithmeticPuzzleEvaluator(dataset_input_dir)
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_combi_calcu(response)
        try:
            match = re.search(r'<<<\[(.*?)\]>>>', output_1)
            solution_1 = match.group(1)
        except:
            solution_1 = ''
        extracted_text_1 = solution_1
        parsed_solution_1 = [elem.strip().strip('"\'') for elem in solution_1.split(',')]
        puzzle = puzzles[i]
        True_false_result_1, message_1 = evaluator.verify_solution(puzzle, parsed_solution_1)
        print('Feedback1: ', message_1)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_combi_calcu(original_response)
        try:
            match = re.search(r'<<<\[(.*?)\]>>>', output_2)
            solution_2 = match.group(1)
        except:
            solution_2 = ''
        extracted_text_2 = solution_2
        parsed_solution_2 = [elem.strip().strip('"\'') for elem in solution_2.split(',')]
        True_false_result_2, message_2 = evaluator.verify_solution(puzzle, parsed_solution_2)
        print('Feedback2: ', message_2)
    elif task_name == 'eight_queens':
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_eight_queens(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_eight_queens(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_eight_queens(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_eight_queens(extracted_text_2, solution_data)
    elif task_name == 'synthesis_decomposition':
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_syn_decom(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_syn_decom(extracted_text_1, solution_data,
                                                                     solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_syn_decom(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_syn_decom(extracted_text_2, solution_data,
                                                                     solution_data['complexity'])
    elif task_name == 'mahjong_pattern':
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_mahjong(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_mahjong(extracted_text_1, solution_data,
                                                                   solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_mahjong(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_mahjong(extracted_text_2, solution_data,
                                                                   solution_data['complexity'])
    elif task_name == 'statistical_counting':
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_stat_counting(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_stat_counting(extracted_text_1, solution_data,
                                                                         solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_stat_counting(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_stat_counting(extracted_text_2, solution_data,
                                                                         solution_data['complexity'])
    elif task_name == 'new_operator':
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_new_op(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_new_op(extracted_text_1, solution_data,
                                                                  solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_new_op(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_new_op(extracted_text_2, solution_data,
                                                                  solution_data['complexity'])
    elif task_name == 'light_puzzles':
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_light(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_light(extracted_text_1, solution_data,
                                                                 solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_light(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_light(extracted_text_2, solution_data,
                                                                 solution_data['complexity'])
    elif task_name == 'reversi':
        # import pdb; pdb.set_trace()
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_reversi(response.replace('\n', ''))
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_reversi(extracted_text_1, solution_data,
                                                                   solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_reversi(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_reversi(extracted_text_2, solution_data,
                                                                   solution_data['complexity'])
    elif task_name == 'matrix_transformation':
        # import pdb; pdb.set_trace()
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_matrix_trans(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_matrix_trans(extracted_text_1, solution_data,
                                                                        solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_matrix_trans(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_matrix_trans(extracted_text_2, solution_data,
                                                                        solution_data['complexity'])
    elif task_name == '2048':
        # import pdb; pdb.set_trace()
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_2048(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_2048(extracted_text_1, solution_data,
                                                                solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_2048(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_2048(extracted_text_2, solution_data,
                                                                solution_data['complexity'])
    elif task_name == 'pooling':
        # import pdb; pdb.set_trace()
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_pooling(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_pooling(extracted_text_1, solution_data,
                                                                   solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_pooling(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_pooling(extracted_text_2, solution_data,
                                                                   solution_data['complexity'])
    elif task_name == 'constrained_linear_arrangement':
        # import pdb; pdb.set_trace()
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_constrained(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_constrained(extracted_text_1, solution_data,
                                                                       solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_constrained(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_constrained(extracted_text_2, solution_data,
                                                                       solution_data['complexity'])
    elif task_name == 'logic_puzzle':
        # import pdb; pdb.set_trace()
        solution_data = solution_data_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_logic_puzzle(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_logic_puzzle(extracted_text_1, solution_data,
                                                                        solution_data['complexity'])

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_logic_puzzle(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_logic_puzzle(extracted_text_2, solution_data,
                                                                        solution_data['complexity'])

    elif task_name == 'pattern_recognition':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_pattern_recognition(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_pattern_recognition(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_pattern_recognition(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_pattern_recognition(extracted_text_2, solution_data)
    elif task_name == 'string_insertion':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_string_insertion(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_string_insertion(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_string_insertion(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_string_insertion(extracted_text_2, solution_data)
    elif task_name == 'letter_logic_diagram':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_letter_logic_diagram(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_letter_logic_diagram(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_letter_logic_diagram(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_letter_logic_diagram(extracted_text_2, solution_data)
    elif task_name == 'string_synthesis':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_string_synthesis(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        # print(f'response: {response}')
        True_false_result_1, message_1 = validate_solution_string_synthesis(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_string_synthesis(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        # print(f'original_response: {original_response}')
        True_false_result_2, message_2 = validate_solution_string_synthesis(extracted_text_2, solution_data)
    elif task_name == 'standard_sudoku':
        question_data = question_matrix_list[i]
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_standard_sudoku(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        # print(f'response: {response}')
        True_false_result_1, message_1 = validate_solution_standard_sudoku(extracted_text_1, question_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_standard_sudoku(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        # print(f'original_response: {original_response}')
        True_false_result_2, message_2 = validate_solution_standard_sudoku(extracted_text_2, question_data)
    elif task_name == 'permutations_and_combinations':
        question_data = question_constrained_list[i]
        # solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_permutations_and_combinations(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        # print(f'response: {response}')
        True_false_result_1, message_1 = validate_solution_permutations_and_combinations(extracted_text_1,
                                                                                         question_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_permutations_and_combinations(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, message_2 = validate_solution_permutations_and_combinations(extracted_text_2,
                                                                                         question_data)
    elif task_name == 'string_deletion_and_modification':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_string_deletion_and_modification(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_string_deletion_and_modification(extracted_text_1,
                                                                                            solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_string_deletion_and_modification(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        # print(f'original_response: {original_response}')
        True_false_result_2, message_2 = validate_solution_string_deletion_and_modification(extracted_text_2,
                                                                                            solution_data)
    elif task_name == 'minesweeper':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_minesweeper(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_minesweeper(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_minesweeper(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        # print(f'original_response: {original_response}')
        True_false_result_2, message_2 = validate_solution_minesweeper(extracted_text_2, solution_data)
    elif task_name == 'cryptanalysis':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_cryptanalysis(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_cryptanalysis(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_cryptanalysis(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        # print(f'original_response: {original_response}')
        True_false_result_2, message_2 = validate_solution_cryptanalysis(extracted_text_2, solution_data)
    elif task_name == 'string_splitting':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_string_splitting(response)
        extracted_text_1, _ = extract_and_check(output_1)
        solution_1 = extracted_text_1
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, message_1 = validate_solution_string_splitting(extracted_text_1, solution_data)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_string_splitting(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        solution_2 = extracted_text_2
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        # print(f'original_response: {original_response}')
        True_false_result_2, message_2 = validate_solution_string_splitting(extracted_text_2, solution_data)
    elif task_name == 'game24':
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_game24(response)
        extracted_text_1, _ = extract_and_check(output_1)
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1 = validate_solution_game24(number_list_item, extracted_text_1)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_game24(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2 = validate_solution_game24(number_list_item, extracted_text_2)
        solution_1 = extracted_text_1; solution_2 = extracted_text_2
    elif task_name == 'letters':
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_letters(response)
        extracted_text_1, _ = extract_and_check(output_1)
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')
        True_false_result_1, explanation_1 = evaluate_response_letters(word, letter, extracted_text_1)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_letters(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')
        True_false_result_2, explanation_2 = evaluate_response_letters(word, letter, extracted_text_2)

        solution_1 = extracted_text_1;
        solution_2 = extracted_text_2
    elif task_name == 'number_multiply':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_number_multiply(response)
        extracted_text_1, _ = extract_and_check(output_1)
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_number_multiply(original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')

        True_false_result_1 = str(solution_data) in extracted_text_1 or str(
            format_number_with_commas(int(solution_data))) in extracted_text_1
        True_false_result_2 = str(solution_data) in extracted_text_2 or str(
            format_number_with_commas(int(solution_data))) in extracted_text_2

        solution_1 = extracted_text_1;
        solution_2 = extracted_text_2
    elif task_name == 'big_bench_hard':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = is_equiv_func_big_bench_hard(question, solution_data, response)
        extracted_text_1, _ = extract_and_check(output_1)
        extracted_text_1 = extracted_text_1.strip()
        print(f'extracted_text_1: {extracted_text_1}')

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = is_equiv_func_big_bench_hard(question, solution_data, original_response)
        extracted_text_2, _ = extract_and_check(output_2)
        extracted_text_2 = extracted_text_2.strip()
        print(f'extracted_text_2: {extracted_text_2}')

        True_false_result_1 = 'Correct' in extracted_text_1 or 'correct' in extracted_text_1
        True_false_result_2 = 'Correct' in extracted_text_2 or 'correct' in extracted_text_2

        solution_1 = extracted_text_1;
        solution_2 = extracted_text_2
    elif task_name == 'gsm':
        solution_data = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_gsm(response)
        extracted_text_1, _ = extract_and_check(output_1)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_gsm(original_response)
        extracted_text_2, _ = extract_and_check(output_2)

        True_false_result_1 = is_equiv_func_gsm(solution_data, extracted_text_1)
        True_false_result_1, _ = extract_and_check(True_false_result_1)
        True_false_result_2 = is_equiv_func_gsm(solution_data, extracted_text_2)
        True_false_result_2, _ = extract_and_check(True_false_result_2)
        solution_1 = extracted_text_1;
        solution_2 = extracted_text_2
    elif task_name == 'math_geometry':
        target_answer = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_math_geometry(response)
        extracted_text_1, _ = extract_and_check(output_1)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_math_geometry(original_response)
        extracted_text_2, _ = extract_and_check(output_2)

        True_false_result_1 = is_equiv_func_math_geometry(target_answer, extracted_text_1)
        True_false_result_1, _ = extract_and_check(True_false_result_1)
        True_false_result_2 = is_equiv_func_math_geometry(target_answer, extracted_text_2)
        True_false_result_2, _ = extract_and_check(True_false_result_2)
        solution_1 = extracted_text_1;
        solution_2 = extracted_text_2
    elif task_name == 'math_counting_and_probability':
        target_answer = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_math_counting_and_probability(response)
        extracted_text_1, _ = extract_and_check(output_1)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_math_counting_and_probability(original_response)
        extracted_text_2, _ = extract_and_check(output_2)

        True_false_result_1 = is_equiv_func_math_counting_and_probability(target_answer, extracted_text_1)
        True_false_result_1, _ = extract_and_check(True_false_result_1)
        True_false_result_2 = is_equiv_func_math_counting_and_probability(target_answer, extracted_text_2)
        True_false_result_2, _ = extract_and_check(True_false_result_2)
        solution_1 = extracted_text_1;
        solution_2 = extracted_text_2
    elif task_name == 'BoxNet1':
        pg_dict = solution_list[i]
        output_1 = None;
        iteration_num_1 = 0
        while output_1 == None and iteration_num_1 < 3:
            iteration_num_1 += 1
            output_1 = extract_equation_with_GPT4_boxnet1(response)
        extracted_text_1, _ = extract_and_check(output_1)

        output_2 = None;
        iteration_num_2 = 0
        while output_2 == None and iteration_num_2 < 3:
            iteration_num_2 += 1
            output_2 = extract_equation_with_GPT4_boxnet1(original_response)
        extracted_text_2, _ = extract_and_check(output_2)

        remaining_box_dict_1, success_failure_1 = score_in_training_set(pg_dict, extracted_text_1)
        remaining_box_dict_2, success_failure_2 = score_in_training_set(pg_dict, extracted_text_2)
        True_false_result_1 = success_failure_1 == 'success'
        True_false_result_2 = success_failure_2 == 'success'

        solution_1 = extracted_text_1;
        solution_2 = extracted_text_2

    print(f'True_false_result from response: {True_false_result_1}')
    print(f'True_false_result from original_response: {True_false_result_2}')
    print(f'target_solution: {solution}')
    print(f'target_answer: {target}')
    print(f'extracted_text from response: {solution_1}')
    print(f'extracted_text from original_response: {solution_2}')
    with open(save_code_dir + f"/True_false_result_1.txt", "w") as f:
        f.write(str(True_false_result_1))
    with open(save_code_dir + f"/True_false_result_2.txt", "w") as f:
        f.write(str(True_false_result_2))
    with open(save_code_dir + f"/extracted_answer_1.txt", "w") as f:
        f.write(extracted_text_1)
    with open(save_code_dir + f"/extracted_answer_2.txt", "w") as f:
        f.write(extracted_text_2)

    if True_false_result_1 == False and True_false_result_2 == False:
        print('False')
        with open(save_code_dir + f"/success_failure.txt", "w") as f:
            f.write('False')
    elif True_false_result_1 == True or True_false_result_2 == True:
        print('True')
        with open(save_code_dir + f"/success_failure.txt", "w") as f:
            f.write('True')
    else:
        raise ValueError('Error: No True_false_result found')

    return True_false_result_1, True_false_result_2


##### Combinatorial Calculation #####
@dataclass
class PuzzleInstance:
    question: str
    numbers: List[int]
    target: int
    solution: List[str]
    complexity: int
    sample_id: int

class ArithmeticPuzzleEvaluator:
    def __init__(self, dataset_dir: str):
        self.dataset_dir = dataset_dir
        
    def read_dataset(self) -> List[PuzzleInstance]:
        """Read all puzzle instances from the dataset directory in sample_i order"""
        puzzles = []
        sample_id = 0
        
        while True:
            sample_dir = os.path.join(self.dataset_dir, f'sample_{sample_id}')
            if not os.path.exists(sample_dir):
                break
                
            try:
                # Read question
                with open(os.path.join(sample_dir, 'question.txt'), 'r') as f:
                    question = f.read().strip()
                
                # Read solution
                with open(os.path.join(sample_dir, 'solution.json'), 'r') as f:
                    solution_data = json.load(f)
                
                puzzle = PuzzleInstance(
                    question=question,
                    numbers=solution_data['numbers'],
                    target=solution_data['target'],
                    solution=solution_data['solution'],
                    complexity=solution_data['complexity'],
                    sample_id=sample_id
                )
                puzzles.append(puzzle)
                sample_id += 1
                
            except Exception as e:
                print(f"Error reading sample_{sample_id}: {e}")
                break
        
        return puzzles

    def evaluate_expression(self, expression: List[str]) -> Optional[float]:
        """Evaluate an arithmetic expression represented as a list of strings"""
        try:
            # Join the expression list into a string
            expr_str = ''.join(expression)
            # Calculate the result
            result = eval(expr_str)
            # Check if result is effectively an integer
            if isinstance(result, (int, float)) and math.isclose(result, round(result), rel_tol=1e-9):
                return round(result)
            return result
        except:
            return None

    def parse_llm_answer(self, answer: str) -> Optional[List[str]]:
        """Parse the LLM's answer from format <<<[symbols]>>> into a list of strings"""
        try:
            # Extract content between <<< and >>>
            match = re.search(r'<<<\[(.*?)\]>>>', answer)
            if not match:
                return None
                
            content = match.group(1)
            # Split by commas and clean up each element
            elements = [elem.strip().strip('"\'') for elem in content.split(',')]
            return elements
            
        except Exception as e:
            print(f"Error parsing LLM answer: {e}")
            return None

    def verify_solution(self, puzzle: PuzzleInstance, parsed_solution) -> Tuple[bool, str]:
        """
        Verify if the LLM's solution is correct for the arithmetic puzzle.

        Args:
            puzzle: PuzzleInstance containing the original problem and solution
            llm_answer: String response from the LLM in format <<<[elements]>>>

        Returns:
            Tuple of (is_correct: bool, message: str)
        """
        # Extract numbers and operations
        numbers = []
        operations = []

        for elem in parsed_solution:
            elem = elem.strip()
            # Try to convert to number
            try:
                num = int(elem)
                numbers.append(num)
            except ValueError:
                # If not a number, must be an operation or parenthesis
                if elem in ['+', '-', '*', '/', '(', ')']:
                    operations.append(elem)
                else:
                    return False, f"Invalid symbol found: {elem}"

        # Check numbers
        if sorted(numbers) != sorted(puzzle.numbers):
            return False, f"Numbers in solution {sorted(numbers)} don't match required numbers {sorted(puzzle.numbers)}"

        # Check operations are valid
        valid_ops = set(['+', '-', '*', '/', '(', ')'])
        if not all(op in valid_ops for op in operations):
            invalid_ops = [op for op in operations if op not in valid_ops]
            return False, f"Invalid operations found: {invalid_ops}"

        # Check parentheses are balanced
        paren_count = 0
        for op in operations:
            if op == '(':
                paren_count += 1
            elif op == ')':
                paren_count -= 1
            if paren_count < 0:  # Closing before opening
                return False, "Unbalanced parentheses - closing before opening"
        if paren_count != 0:  # Unequal number of open/close
            return False, "Unbalanced parentheses - unclosed parentheses"

        # Evaluate the expression
        try:
            result = self.evaluate_expression(parsed_solution)
            if result is None:
                return False, "Could not evaluate expression"

            if not math.isclose(result, puzzle.target, rel_tol=1e-9):
                return False, f"Expression evaluates to {result}, target is {puzzle.target}"

        except Exception as e:
            return False, f"Error evaluating expression: {str(e)}"

        # If we got here, everything checked out
        return True, "Correct solution"

def extract_equation_with_GPT4_combi_calcu(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<list of integers and symbols>>>, like <<<[(, 6, /, 1, ), +, 4]>>>, or <<<[6, +, 3, +, 1]>>>, or <<<[9, *, (, 6, -, 9, /, 4, )]>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer empty list <<<[0]>>. Symbols and numbers cannot be connected into one item, items like (6, 3), (4 are totally wrong. If finded, correct it like correcting wrong [(6, +, 3)] to correct [(, 6, +, 3, )].\n' \
             'Also complete expression like (9 + 9) + (4 + 6) should be corrected to list of items [(, 9, +, 9, ), +, (, 4, +, 6, )]\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

##### Logical Equation #####
def read_dataset_logical_equation(dataset_input_dir):
    question_list = []; solution_list = []
    for num_letters, num_constraints, values in [
        (9, 7, [1, 3, 4, 9, 16, 27, 36, 80, 121]),
        (9, 8, [3, 6, 9, 20, 32, 36, 80, 121, 120]),
        (11, 10, [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]),
        (11, 11, [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]),
        (13, 10, [1, 2, 3, 5, 7, 16, 15, 24, 10, 45, 28, 36, 50]),
        (13, 11, [2, 3, 5, 7, 16, 15, 24, 10, 45, 28, 36, 50, 96]),
        (13, 12, [2, 3, 5, 7, 16, 15, 24, 10, 45, 28, 36, 50, 96]),
        (15, 12, [2, 3, 5, 7, 16, 15, 24, 10, 45, 28, 36, 50, 78, 90, 100]),
        (15, 13, [1, 3, 5, 7, 16, 15, 24, 12, 45, 34, 36, 56, 78, 95, 100]),
        (15, 14, [1, 3, 5, 7, 16, 15, 24, 12, 45, 34, 36, 56, 78, 95, 100]),
    ]:
        for i in range(15):
            base_dir = os.path.join(dataset_input_dir, f'sample_{num_letters}_{num_constraints}_{i}')
            with open(os.path.join(base_dir, f'question.txt'), 'r') as f:
                question = f.read()
                question_list.append(question)
            with open(os.path.join(base_dir, f'solution.json'), 'r') as f:
                solution = json.load(f)
                solution_list.append(solution)

    # Create list of indices
    indices = list(range(len(solution_list)))
    # Shuffle the indices
    #random.shuffle(indices)

    # Create new lists using the shuffled indices
    shuffled_solutions = [solution_list[i] for i in indices]
    shuffled_questions = [question_list[i] for i in indices]
    return shuffled_solutions, shuffled_questions

def extract_equation_with_GPT4_logical_equation(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<list of integers>>>, like <<<[1, 4, 20, 6]>>>, or <<<[20, 43, 12, 50, 64]>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer empty list <<<[]>>.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def verify_solution_logical_equation(example_solution: List[int], proposed_solution: List[int], constraints: List[str]) -> bool:
    """Verify if the proposed solution satisfies all constraints"""

    if len(proposed_solution) != len(example_solution):
        return False
    if None in proposed_solution:
        return False
    #if all(isinstance(item, int) for item in proposed_solution):
    #    return False
    letters = [chr(65 + i) for i in range(len(proposed_solution))]
    solution_dict = {letters[i]: proposed_solution[i] for i in range(len(proposed_solution))}

    print(f'\nexample_solution: {example_solution}')
    print(f'proposed_solution: {proposed_solution}\n')
    for constraint in constraints:
        print(constraint)
    print('\n')

    def parse_term(term: str) -> float:
        if term[0].isalpha():  # Single letter
            return solution_dict[term]
        else:  # Term like '3B'
            multiplier = float(term[:-1])
            letter = term[-1]
            return multiplier * solution_dict[letter]

    for constraint in constraints:
        parts = constraint.split()

        if len(parts) == 3:  # A = B, A > B, A < B
            left_val = parse_term(parts[0])
            right_val = parse_term(parts[2])

            if parts[1] == "=":
                if left_val != right_val:
                    return False
            elif parts[1] == ">":
                if left_val <= right_val:
                    return False
            elif parts[1] == "<":
                if left_val >= right_val:
                    return False

        elif len(parts) == 5:  # A + B = 5, B - D = 1
            left_val = parse_term(parts[0])
            right_val = parse_term(parts[2])
            result = float(parts[4])

            if parts[1] == "+":
                if left_val + right_val != result:
                    return False
            elif parts[1] == "-":
                if left_val - right_val != result:
                    return False
    for item in proposed_solution:
        if int(item) not in example_solution:
            #print(f"Item {item} from proposed not in example")
            return False
    # Check if all numbers in proposed solution are valid
    #return sorted(proposed_solution) == sorted(example_solution)
    return True

##### Eight Queens #####
def validate_solution_eight_queens(response: str, solution_data: Dict) -> Tuple[bool, str]:
        # Parse response into queen positions
        positions = response.split(',')
        queens = []
        try:
            for pos in positions:
                row, col = map(int, pos.strip().split())
                queens.append((row, col))
        except:
            pass

        # Check number of queens
        if len(queens) != 8:
            return False, f"Wrong number of queens: {len(queens)}/8"

        # Check for blocked positions
        blocked = solution_data['blocked_positions']
        for queen in queens:
            if list(queen) in blocked:
                return False, f"Queen placed on blocked position: {queen}"

        # Check for conflicts
        for i, q1 in enumerate(queens):
            for q2 in queens[i + 1:]:
                if (q1[0] == q2[0] or  # Same row
                        q1[1] == q2[1] or  # Same column
                        abs(q1[0] - q2[0]) == abs(q1[1] - q2[1])):  # Same diagonal
                    return False, f"Conflicting queens: {q1} and {q2}"

        return True, "Valid solution"


def read_dataset_eight_queens(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_eight_queens(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<position pairs of integers>>>, like <<<0 3, 1 0, 2 4>>>, or <<<2 0, 4 3, 1 2, 5 0, 6 4>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0 0>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation


##### Synthesis and decomposition #####
def validate_solution_syn_decom(response: str, solution: Dict[str, Union[List[str], List[List[str]]]],
                                complexity: int) -> Tuple[bool, str]:
    # import pdb; pdb.set_trace()
    # pattern = r'<<<\s*\[(.*?)\]\s*>>>'
    # match = re.search(pattern, response, re.DOTALL)
    # if not match:
    #     return False, f"Wrong response fromat"

    # answer_str = match.group(1)
    # Split by comma and clean up each element
    answer_parts = [part.strip() for part in response.split(',')]
    # import pdb; pdb.set_trace()

    # # Check required fields
    # if 'answer' not in llm_solution or 'process' not in llm_solution:
    #     return False, f"JSON needs to include answer and process fields"

    # # Check answer format
    # if not isinstance(llm_solution['answer'], list):
    #     return False, f"answer field should contain a list"

    # # Check process format
    # if not isinstance(llm_solution['process'], list):
    #     return False, f"process field should contain a list"

    # Check number of elements based on complexity
    expected_length = 5 if complexity > 4 else (4 if complexity > 2 else 3)
    if len(answer_parts) != expected_length:
        return False, f"number of agricultural product does not match the question"

    # Check if answer matches
    if complexity > 2:
        if answer_parts != solution['solution']['answer'][:expected_length]:
            return False, f"answer is not correct"
    else:
        if answer_parts[:2] != solution['solution']['answer'][:2] or answer_parts[-1] != solution['solution']['answer'][
            -2]:
            return False, f"answer is not correct"

    # # Check if process steps are valid
    # for step in llm_solution['process']:
    #     if not isinstance(step, list) or len(step) != expected_length:
    #         return False, f"process is not valid"

    return True, "Valid solution"

def read_dataset_syn_decom(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_syn_decom(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<list of values of crops>>>, like <<<3, 1, 2>>>, or <<<1, 3, 1, 2, 0>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation


##### Mahjong-type #####
def validate_solution_mahjong(response: str, solution: Dict[str, Union[List[str], List[List[str]]]], complexity: int) -> \
Tuple[bool, str]:
    # import pdb; pdb.set_trace()
    # answer_parts = [part.strip() for part in response.split(',')]
    if response == '':
        return False, f"answer cannot be empty"

    elif int(response) != solution['solution']:
        return False, f"answer is not correct"

    return True, "Valid solution"


def read_dataset_mahjong(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_mahjong(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<number in final round>>>, like <<<2>>>, or <<<0>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<-1>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### Statistical counting #####
def validate_solution_stat_counting(response: str, solution: Dict[str, Union[List[str], List[List[str]]]],
                                    complexity: int) -> Tuple[bool, str]:
    # import pdb; pdb.set_trace()
    # answer_parts = [part.strip() for part in response.split(',')]
    if response == '':
        return False, f"answer cannot be empty"

    elif int(response) != solution['solution']:
        return False, f"answer is not correct"

    return True, "Valid solution"


def read_dataset_stat_counting(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_stat_counting(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<final score>>>, like <<<2>>>, or <<<5>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### new operator #####
def validate_solution_new_op(response: str, solution: Dict[str, Union[List[str], List[List[str]]]], complexity: int) -> \
Tuple[bool, str]:
    # import pdb; pdb.set_trace()
    # answer_parts = [part.strip() for part in response.split(',')]
    try:
        if response == '':
            return False, f"answer cannot be empty"

        elif float(response) != float(solution['solution']):
            return False, f"answer is not correct"
    except:
        return False, f"answer cannot be convert to float"

    return True, "Valid solution"


def read_dataset_new_op(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_new_op(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<final result>>>, like <<<2>>>, or <<<5>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### light out #####
def validate_solution_light(response: str, solution: Dict[str, Union[List[str], List[List[str]]]], complexity: int) -> \
Tuple[bool, str]:
    try:
        answer_parts = [int(part.strip()) for part in response.split(',')]
        # import pdb; pdb.set_trace()

        if len(answer_parts) != solution['n'] ** 2:
            return False, f"number of lights in light network is not correct"

        # Check if answer matches
        if answer_parts != solution['solution']:
            return False, f"answer is not correct"
    except:
        return False, f"answer format is not correct"

    return True, "Valid solution"


def read_dataset_light(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_light(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<final light network>>>, like <<<1, 1, 0, 0>>>, or <<<1, 0, 1, 0, 1, 1, 1, 1, 0>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### reversi #####
def validate_solution_reversi(response: str, solution: Dict[str, Union[List[str], List[List[str]]]], complexity: int) -> \
Tuple[bool, str]:
    try:
        # import pdb; pdb.set_trace()
        answer_parts = [part.strip() for part in response.split(',')]
        solution_answer = [part.strip() for part in solution['solution'].split(',')]

        if len(answer_parts) != solution['board_size'] ** 2:
            return False, f"number of position in the board is not correct"

        # Check if answer matches
        if answer_parts != solution_answer:
            return False, f"answer is not correct"
    except:
        return False, f"answer format is not correct"

    return True, "Valid solution"


def read_dataset_reversi(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_reversi(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<final board situations>>>, like <<<*, *, *, *, 1, 1, 0, 0, *, *, *, *, 1, 1, 0, 0>>>, or <<<1, 0, 1, 0, 1, 1, 1, 1, 0>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Ignore any newline character in the answer. Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### matrix transformation #####
def validate_solution_matrix_trans(response: str, solution: Dict[str, Union[List[str], List[List[str]]]],
                                   complexity: int) -> Tuple[bool, str]:
    try:
        answer_parts = [part.strip() for part in response.split(',')]
        solution = [part.strip() for part in solution['formatted_solution'].split(',')]
        # import pdb; pdb.set_trace()

        if len(answer_parts) != len(solution):
            return False, f"number of position in the board is not correct"

        # Check if answer matches
        if answer_parts != solution:
            return False, f"answer is not correct"
    except:
        return False, f"answer format is not correct"

    return True, "Valid solution"


def read_dataset_matrix_trans(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_matrix_trans(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<final board situations>>>, like <<<A, B, C, D>>>, or <<<A, B, C, D, E, F, H, I, G>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### 2048 #####
def validate_solution_2048(response: str, solution: Dict[str, Union[List[str], List[List[str]]]], complexity: int) -> \
Tuple[bool, str]:
    try:
        # import pdb; pdb.set_trace()
        answer_parts = [part.strip() for part in response.split(',')]
        solution = [part.strip() for part in solution['solution'].replace('\n', ',').split(',')]
        # import pdb; pdb.set_trace()

        if len(answer_parts) != len(solution):
            return False, f"number of position in the board is not correct"

        # Check if answer matches
        if answer_parts != solution:
            return False, f"answer is not correct"
    except:
        return False, f"answer format is not correct"

    return True, "Valid solution"


def read_dataset_2048(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_2048(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<final board layout>>>, like <<<2, 2, 0, 8>>>, or <<<2, 4, 8, 16, 2, 4, 0, 0, 0>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### pooling #####
def validate_solution_pooling(response: str, solution: Dict[str, Union[List[str], List[List[str]]]], complexity: int) -> \
Tuple[bool, str]:
    try:
        # import pdb; pdb.set_trace()
        answer_parts = [int(part.strip()) for part in response.split(',')]
        solution = np.concatenate(solution['solution']).tolist()
        # import pdb; pdb.set_trace()

        if len(answer_parts) != len(solution):
            return False, f"number of position in the board is not correct"

        # Check if answer matches
        if answer_parts != solution:
            return False, f"answer is not correct"
    except:
        return False, f"answer format is not correct"

    return True, "Valid solution"


def read_dataset_pooling(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_pooling(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<pooling result>>>, like <<<2, 2, 0, 8>>>, or <<<2, 4, 8, 16, 2, 4, 0, 0, 0>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                            user_prompt_list=[prompt + response], response_total_list=[],
                                            logprobs=False)
    return extract_equation


##### Constrained linear arrangement #####
def validate_solution_constrained(response: str, solution: Dict[str, Union[List[str], List[List[str]]]],
                                  complexity: int) -> Tuple[bool, str]:
    try:
        # import pdb; pdb.set_trace()
        answer_parts = [part.strip() for part in response.split(',')]
        solution_answer = solution['solution']  # [part.strip() for part in solution['solution'].split(',')]

        if len(answer_parts) != len(solution_answer):
            return False, f"number of played cards needs to be same as number of rounds"

        elements = {
            "pieces": ["A", "B", "C", "D", "E"],  # Metal, Wood, Water, Fire, Earth
            "generates": {"A": "B", "B": "C", "C": "D", "D": "E", "E": "A"},
            "overcomes": {"A": "C", "C": "E", "E": "B", "B": "D", "D": "A"}
        }
        stick = {
            "pieces": ["A", "B", "C", "D"],  # Stick, Tiger, Chicken, Worm
            "beats": {"A": "B", "B": "C", "C": "D", "D": "A"}
        }
        animal = {
            "pieces": ["A", "B", "C", "D", "E", "F", "G", "H"],  # Elephant to Mouse
            "special": {"H": "A"},  # Mouse beats Elephant
            "hierarchy": {"A": 8, "B": 7, "C": 6, "D": 5, "E": 4, "F": 3, "G": 2, "H": 1}
        }
        if solution['game_type'] == "elements":
            for i, result in enumerate(solution['results']):
                player_move = solution['player_moves'][i]
                if result == 'draw':
                    if answer_parts[i] == elements['generates'][player_move] or answer_parts[i] == \
                            elements['overcomes'][player_move]:
                        return False, f"answered card cannot result in a draw for round {i + 1}"
                if result == 'win':
                    if not answer_parts[i] == solution_answer[i]:
                        return False, f"answered card cannot result in a win for round {i + 1}"
                if result == 'loss':
                    if not answer_parts[i] == solution_answer[i]:
                        return False, f"answered card cannot result in a loss for round {i + 1}"
        elif solution['game_type'] == "stick_game":
            for i, result in enumerate(solution['results']):
                player_move = solution['player_moves'][i]
                if result == 'draw':
                    if answer_parts[i] == stick['beats'][player_move] or player_move == stick['beats'][answer_parts[i]]:
                        return False, f"answered card cannot result in a draw for round {i + 1}"
                if result == 'win':
                    if not answer_parts[i] == solution_answer[i]:
                        return False, f"answered card cannot result in a win for round {i + 1}"
                if result == 'loss':
                    if not answer_parts[i] == solution_answer[i]:
                        return False, f"answered card cannot result in a loss for round {i + 1}"
        elif solution['game_type'] == "animal":
            if len(answer_parts) != len(set(answer_parts)):
                return False, f"answered card cannot repeat"
            for i, result in enumerate(solution['results']):
                player_move = solution['player_moves'][i]
                if result == 'draw':
                    if not answer_parts[i] == player_move:
                        return False, f"answered card cannot result in a draw for round {i + 1}"
                if result == 'win':
                    if not (animal["hierarchy"][player_move] > animal["hierarchy"][answer_parts[i]] or (
                            player_move == "H" and answer_parts[i] == "A")):
                        return False, f"answered card cannot result in a win for round {i + 1}"
                if result == 'loss':
                    if not (animal["hierarchy"][player_move] < animal["hierarchy"][answer_parts[i]] or (
                            player_move == "A" and answer_parts[i] == "H")):
                        return False, f"answered card cannot result in a loss for round {i + 1}"
    except:
        return False, f"answer format is not correct"

    return True, "Valid solution"


def read_dataset_constrained(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break
    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_constrained(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<opponent played card result>>>, like <<<A, B, C>>>, or <<<C, D, E, F>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                    user_prompt_list=[prompt + response], response_total_list=[],
                                    logprobs=False)
    return extract_equation


##### Logic puzzle  #####
def validate_solution_logic_puzzle(response: str, solution_data: Dict, complexity: int) -> Tuple[bool, str]:
    # Parse response into queen positions
    positions = response.split(',')
    numbers = []
    try:
        for pos in positions:
            row, col = map(int, pos.strip().split())
            numbers.append((row, col))
    except:
        pass
    # import pdb; pdb.set_trace()
    solution = solution_data['solution']
    grid = solution_data['grid']
    # Check number of queens
    if len(numbers) != len(solution):
        return False, f"Wrong number of selected number positions"

    for row in range(len(grid)):
        summation = 0
        product = 1
        for column in range(len(grid)):
            if (row, column) in numbers:
                summation += grid[row][column]
                product *= grid[row][column]
        if solution_data['constraints']['type'] == 'sum_le_4':
            if summation > 4:
                return False, f"the constraint of sum of the selected numbers in each row should be less than or equal to 4 is not fulfilled"
        if solution_data['constraints']['type'] == 'product_gt_0':
            if product <= 0:
                return False, f"the constraint of product of the selected numbers in each row should be greater than 0 is not fulfilled"

    for column in range(len(grid)):
        summation = 0
        product = 1
        for row in range(len(grid)):
            if (row, column) in numbers:
                summation += grid[row][column]
                product *= grid[row][column]
        if solution_data['constraints']['type'] == 'sum_le_4':
            if summation > 4:
                return False, f"the constraint of sum of the selected numbers in each column should be less than or equal to 4 is not fulfilled"
        if solution_data['constraints']['type'] == 'product_gt_0':
            if product <= 0:
                return False, f"the constraint of product of the selected numbers in each column should be greater than 0 is not fulfilled"

    return True, "Valid solution"


def read_dataset_logic_puzzle(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break
    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_logic_puzzle(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<selected position pairs>>>, like <<<0 3, 1 0, 2 4>>>, or <<<2 0, 4 3, 1 2, 5 0, 6 4>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0 0>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: '
    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False,
                                    user_prompt_list=[prompt + response], response_total_list=[],
                                    logprobs=False)
    return extract_equation


#####Pattern Recognition#######
def read_dataset_pattern_recognition(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles


def extract_equation_with_GPT4_pattern_recognition(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<one position pairs of integers (row, column)>>>, like <<<1,3>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<0,0>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def validate_solution_pattern_recognition(response: str, solution_data: list) -> Tuple[bool, str]:
    try:
        # Parse response into position
        pos = response.split(',')
        row, col = pos[0].strip(), pos[1].strip()
        print(row, col)

        row1, col1 = solution_data[0], solution_data[1]
        print(row1, col1)

        # Check for conflicts
        if not(int(row) == int(row1) and int(col) == int(col1)):
            return False, f"Conflicting position: {(row, col)} and {(row1, col1)}"

        return True, "Valid solution"
    except:
        return False, f"answer format is not correct"

#####String Insertion#######
def read_dataset_string_insertion(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_string_insertion(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<a string>>>, like <<<ABCDEACE>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def validate_solution_string_insertion(response: str, solution_data: str) -> Tuple[bool, str]:
    try:
        # Parse response into position
        str1 = str(response.strip())
        print(str1)

        str2 = solution_data
        print(str2)

        # Check for conflicts
        if str1 != str2:
            return False, f"Conflicting string: {str1} and {str2}"

        return True, "Valid solution"
    except:
        return False, f"answer format is not correct"


#####Letter Logic Diagram#######
def read_dataset_letter_logic_diagram(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break
    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_letter_logic_diagram(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<a 7*7 matrix>>>, like <<<[[a,b,c,d,e,f,g],[b,.......], ......]>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<[]>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def validate_solution_letter_logic_diagram(
        response: str,
        solution_data: List[List[str]]
) -> Tuple[bool, str]:
    try:
        # 1. Safely parse the response string into a Python object
        inner_str = response.strip('[]')
        sub_strs = inner_str.split('],[')

        parsed_response = []
        for part in sub_strs:
            clean_part = part.strip('[]')
            items = clean_part.split(',')
            items = [x.strip() for x in items]
            parsed_response.append(items)
        print(parsed_response)
        print(solution_data)
        # 2. Validate the data structure
        # Check if parsed_response is a list of length 7
        if len(parsed_response) != 7:
            return False, "The response data is not a list of length 7 (rows)."

        # For each row, check if it is a list of length 7 and contains single-character strings
        for i, row in enumerate(parsed_response):
            if len(row) != 7:
                return False, f"Row {i + 1} is not a list of length 7."
            for j, item in enumerate(row):
                if parsed_response[i][j] != solution_data[i][j]:
                    return False, "The answer does not match the correct solution."

        return True, "Correct answer."
    except:
        return False, f"answer format is not correct"


#####String Synthesis#######
def read_dataset_string_synthesis(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_string_synthesis(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a string of the number of certain block>>>, like <<<101010000>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def validate_solution_string_synthesis(
        response: str,
        solution_data: List[int]
) -> Tuple[bool, str]:

    try:
        print("response_list"+response)
        if response == "":
            return False, f"Empty list"
        response = response.replace("[", "").replace("]", "").replace("'", "").replace(",", "").replace(" ", "")

        try:
            response_list = [int(x) for x in response]
        except ValueError as e:
            return False, f"Conversion error: {e}"

        print(response_list)
        print(solution_data)

        if response_list == solution_data:
            return True, "Valid solution"
        else:
            return False, f"Conflicting list: {response_list} and {solution_data}"
    except:
        return False, f"answer format is not correct"

#####Standard Sudoku#######
def read_dataset_standard_sudoku(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_standard_sudoku(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a 9*9 matrix >>>, like <<<[[1,2,3,4,5,6,7,8,9],[2,...],...,[]]>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def check_sudoku_solution(
        puzzle: List[List[int]],
        solution_candidate: List[List[int]]
) -> bool:
    # Check givens
    for r in range(9):
        for c in range(9):
            if puzzle[r][c] != 0 and puzzle[r][c] != solution_candidate[r][c]:
                return False
    # Check rows
    for r in range(9):
        row_vals = set()
        for c in range(9):
            val = solution_candidate[r][c]
            if val < 1 or val > 9 or val in row_vals:
                return False
            row_vals.add(val)
    # Check columns
    for c in range(9):
        col_vals = set()
        for r in range(9):
            val = solution_candidate[r][c]
            if val < 1 or val > 9 or val in col_vals:
                return False
            col_vals.add(val)
    # Check 3x3 subgrids
    for sub_row in range(0, 9, 3):
        for sub_col in range(0, 9, 3):
            box_vals = set()
            for r in range(sub_row, sub_row + 3):
                for c in range(sub_col, sub_col + 3):
                    val = solution_candidate[r][c]
                    if val in box_vals:
                        return False
                    box_vals.add(val)
    return True

def  validate_solution_standard_sudoku(
        response: str,
        question_data: List[List[int]]
) -> Tuple[bool, str]:

    try:
        if response == "":
            return False, "False answer."
        try:
            parsed_response = json.loads(response)
            is_correct = check_sudoku_solution(question_data, parsed_response)
        except:
            return False, "False answer."

        if is_correct:
            return True, "Correct answer."
        else:
            return False, "False answer."
    except:
        return False, f"answer format is not correct"

#####permutations_and_combinations#######
def read_dataset_permutations_and_combinations(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_permutations_and_combinations(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a list of strings >>>, like <<<["A", "B", "C", "D", "E"]>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<[]>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def parse_constraints(constraints_text: str) -> List[Tuple]:
    """
    Parse each line of constraints_text into a structured list of constraints.
    We handle 5 types of constraints, matching how they're generated in the puzzle code:

      1) ("fixed_position", item, pos)
         e.g. "Book A must be placed in position 3."

      2) ("left_of", item1, item2)
         e.g. "Book A must be to the left of book B."

      3) ("not_in_position", item, pos)
         e.g. "Book C cannot be placed in position 2."

      4) ("right_of", item1, item2)
         e.g. "Book D must be to the right of book A."

      5) ("adjacent_to", item1, item2)
         e.g. "Book E must be adjacent to book F."

    Returns: a list of tuples like [("left_of", "A", "B"), ("fixed_position", "C", 4), ...]
    """

    # Split the text by newlines
    lines = constraints_text.strip().split('\n')
    parsed_constraints = []

    # Regex patterns for each known constraint type
    # 1) fixed_position
    #    "Book X must be placed in position Y."
    pattern_fixed_position = re.compile(
        r'Book\s+([A-Za-z])\s+must\s+be\s+placed\s+in\s+position\s+(\d+)'
    )

    # 2) left_of
    #    "Book X must be to the left of book Y."
    pattern_left_of = re.compile(
        r'Book\s+([A-Za-z])\s+must\s+be\s+to\s+the\s+left\s+of\s+book\s+([A-Za-z])'
    )

    # 3) not_in_position
    #    "Book X cannot be placed in position Y."
    pattern_not_in_pos = re.compile(
        r'Book\s+([A-Za-z])\s+cannot\s+be\s+placed\s+in\s+position\s+(\d+)'
    )

    # 4) right_of
    #    "Book X must be to the right of book Y."
    pattern_right_of = re.compile(
        r'Book\s+([A-Za-z])\s+must\s+be\s+to\s+the\s+right\s+of\s+book\s+([A-Za-z])'
    )

    # 5) adjacent_to
    #    "Book X must be adjacent to book Y."
    pattern_adjacent_to = re.compile(
        r'Book\s+([A-Za-z])\s+must\s+be\s+adjacent\s+to\s+book\s+([A-Za-z])'
    )

    for line in lines:
        line = line.strip()

        # Try each pattern in turn:
        m_fp = pattern_fixed_position.search(line)
        if m_fp:
            item = m_fp.group(1)
            pos = int(m_fp.group(2))
            parsed_constraints.append(("fixed_position", item, pos))
            continue

        m_lo = pattern_left_of.search(line)
        if m_lo:
            item1 = m_lo.group(1)
            item2 = m_lo.group(2)
            parsed_constraints.append(("left_of", item1, item2))
            continue

        m_nip = pattern_not_in_pos.search(line)
        if m_nip:
            item = m_nip.group(1)
            pos = int(m_nip.group(2))
            parsed_constraints.append(("not_in_position", item, pos))
            continue

        m_ro = pattern_right_of.search(line)
        if m_ro:
            item1 = m_ro.group(1)
            item2 = m_ro.group(2)
            parsed_constraints.append(("right_of", item1, item2))
            continue

        m_adj = pattern_adjacent_to.search(line)
        if m_adj:
            item1 = m_adj.group(1)
            item2 = m_adj.group(2)
            parsed_constraints.append(("adjacent_to", item1, item2))
            continue

        # If no known pattern matched:
        print(f"[Warning] Unrecognized constraint format:\n  {line}")

    return parsed_constraints


def check_constraints(arrangement: List[str], constraints: List[Tuple]) -> Tuple[bool, str]:
    """
    Given a final arrangement (e.g. ["A","B","E","C","D"])
    and a list of constraint tuples, verify each constraint.

    Return:
      (True, "OK") if all constraints are satisfied
      (False, <reason>) if any constraint is violated
    """

    # Build a map { item: position_index }, using 1-based indexing
    position_map = {book: i+1 for i, book in enumerate(arrangement)}

    for constraint in constraints:
        ctype = constraint[0]

        if ctype == "fixed_position":
            # e.g. ("fixed_position", "E", 3)
            _, item, pos = constraint
            actual = position_map.get(item, -1)
            if actual != pos:
                return False, (
                    f"Constraint violated: Book {item} must be placed in position {pos}, "
                    f"but is at position {actual}."
                )

        elif ctype == "left_of":
            # e.g. ("left_of", "A", "C")
            _, item1, item2 = constraint
            if position_map.get(item1, 99999) >= position_map.get(item2, -99999):
                return False, (
                    f"Constraint violated: Book {item1} must be to the left of book {item2}, "
                    f"but positions are {item1}={position_map[item1]}, {item2}={position_map[item2]}."
                )

        elif ctype == "not_in_position":
            # e.g. ("not_in_position", "C", 1)
            _, item, pos = constraint
            if position_map.get(item, -1) == pos:
                return False, (
                    f"Constraint violated: Book {item} cannot be placed in position {pos}, "
                    f"but it is at position {pos}."
                )

        elif ctype == "right_of":
            # e.g. ("right_of", "D", "A")
            _, item1, item2 = constraint
            if position_map.get(item1, -1) <= position_map.get(item2, 99999):
                return False, (
                    f"Constraint violated: Book {item1} must be to the right of book {item2}, "
                    f"but positions are {item1}={position_map[item1]}, {item2}={position_map[item2]}."
                )

        elif ctype == "adjacent_to":
            # e.g. ("adjacent_to", "A", "B")
            _, item1, item2 = constraint
            pos_diff = abs(position_map.get(item1, -1) - position_map.get(item2, -1))
            if pos_diff != 1:
                return False, (
                    f"Constraint violated: Book {item1} must be adjacent to book {item2}, "
                    f"positions are {item1}={position_map[item1]}, {item2}={position_map[item2]}."
                )

        else:
            return False, f"Unknown constraint type: {constraint}"

    # If no constraint was violated
    return True, "OK"


def check_solution(constraints_text: str, correct_solution: List[str]) -> Tuple[bool, str]:
    """
    Entry point to validate a solution against puzzle constraints described in text.
    1) parse the constraints
    2) check them one by one on the provided correct_solution

    Returns (True, "OK") if all constraints are satisfied,
    otherwise (False, error_reason).
    """
    constraints = parse_constraints(constraints_text)
    valid, reason = check_constraints(correct_solution, constraints)
    return valid, reason
def  validate_solution_permutations_and_combinations(
        response: str,
        question_data: str
) -> Tuple[bool, str]:

    try:
        print("response" + response)
        print("question_data" + question_data)
        if response == "[]" or response == "":
            return False, "Empty answer"

        response_solution = json.loads(response)
        result, message = check_solution(question_data, response_solution)

        if result:
            return True, message
        else:
            return False, message
    except:
        return False, f"answer format is not correct"


#####String Deletion And Modification#######
def read_dataset_string_deletion_and_modification(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_string_deletion_and_modification(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a string>>>, like <<bbbababab>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def  validate_solution_string_deletion_and_modification(
        response: str,
        solution_data: str
) -> Tuple[bool, str]:
    try:
        print("response" + response)
        print("solution_data" + solution_data)
        if response == "[]" or response == "":
            return False, "Empty answer"

        if str(response) == str(solution_data):
            return True, "Correct Answer"
        else:
            return False, "False"
    except:
        return False, f"answer format is not correct"


#####Minesweeper#######
def read_dataset_minesweeper(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_minesweeper(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a matrix >>>, like <<<[[1,2], [1,3], [2,1], ...]>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<[]>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def  validate_solution_minesweeper(
        response: str,
        solution_data: List[List[int]]
) -> Tuple[bool, str]:
    print("response" + response)
    print("solution_data")
    print(solution_data)

    try:
        if response == "[]" or response == "":
            return False, "Empty answer"

        try:
            list_data = json.loads(response)
        except json.JSONDecodeError as e:
            return False, "JSON Error"

        print(list_data)

        list_set = set(map(tuple, list_data))  # 将每个坐标转换为元组，再转换为集合
        solution_set = set(map(tuple, solution_data))

        if list_set == solution_set:
            print("The coordinates are the same.")
            return True, "Correct Answer"
        else:
            print("The coordinates are different.")
            return False, "False"
    except:
        return False, f"answer format is not correct"

#####String Deletion And Modification#######
def read_dataset_string_deletion_and_modification(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_string_deletion_and_modification(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a string>>>, like <<bbbababab>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def  validate_solution_string_deletion_and_modification(
        response: str,
        solution_data: str
) -> Tuple[bool, str]:
    try:
        print("response" + response)
        print("solution_data" + solution_data)
        if response == "[]" or response == "":
            return False, "Empty answer"

        if str(response) == str(solution_data):
            return True, "Correct Answer"
        else:
            return False, "False"
    except:
        return False, f"answer format is not correct"


#####Cryptanalysis#######
def read_dataset_cryptanalysis(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_cryptanalysis(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a string >>>, like <<<AC10>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def  validate_solution_cryptanalysis(
        response: str,
        solution_data: List[str]
) -> Tuple[bool, str]:
    print("response" + response)
    print("solution_data")
    print(solution_data)

    try:
        if response == "[]" or response == "":
            return False, "Empty answer"
        solution_str = ''.join(solution_data)
        print(solution_str)

        if response == solution_str:
            print("The strings are the same.")
            return True, "Correct Answer"
        else:
            print("The strings are different.")
            return False, "False"
    except:
        return False, f"answer format is not correct"


#####String Splitting#######
def read_dataset_string_splitting(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    while True:
        sample_dir = os.path.join(dataset_dir, f'sample_{sample_id}')
        if not os.path.exists(sample_dir):
            print(f'No sample_{sample_dir} found')
            break

        try:
            # Read question and solution
            question_path = os.path.join(sample_dir, 'question.txt')
            solution_path = os.path.join(sample_dir, 'solution.json')

            with open(question_path, 'r') as f:
                question = f.read().strip()

            with open(solution_path, 'r') as f:
                solution_data = json.load(f)

            puzzles.append({
                'sample_id': sample_id,
                'question': question,
                'solution_data': solution_data
            })
            sample_id += 1

        except Exception as e:
            print(f"Error reading sample_{sample_id}: {e}")
            break

    shuffled_puzzles = puzzles.copy()  # Create a copy to avoid modifying original
    #random.shuffle(shuffled_puzzles)
    return shuffled_puzzles

def extract_equation_with_GPT4_string_splitting(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<< a string representing the outcome in the order of machines A, B, C, ... parts X, Y, Z, ...>>>, like <<<112...>>>,.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<>>>. If the initial answer follows wrong format, then correct it.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def  validate_solution_string_splitting(
        response: str,
        solution_data: List[str]
) -> Tuple[bool, str]:
    print("response" + response)
    print("solution_data")
    print(solution_data)

    try:
        if response == "[]" or response == "":
            return False, "Empty answer"
        solution_str = ''.join(solution_data)
        print(solution_str)

        if response == solution_str:
            print("The strings are the same.")
            return True, "Correct Answer"
        else:
            print("The strings are different.")
            return False, "False"
    except:
        return False, f"answer format is not correct"

##### Game24 #####
def extract_equation_with_GPT4_game24(response):
    prompt = 'Your task is to extract the equation from the given answer by another LLM:\n' \
             'Note that the equation should include four numbers be in the form like ((11 * 8) + 8) / 4, ((3 * 5) - 12) * 8, ((7 - 4) * 11) - 9 = 24,' \
             '((6 / 3) * 7) + 10, (37 - (29 - 16)), ((19 + 18) - 13) = 24\n No other symbols or texts should be included. If included, remove it.' \
             'Here is the reponse, return your answer with the format <<<equation>>>, like <<<((7 - 4) * 11) - 9 = 24>>>. ' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list=[prompt + response], response_total_list=[], logprobs=False)
    return extract_equation

def validate_solution_game24(number_list, extracted_text):
    if "=" not in extracted_text:
        extracted_text = extracted_text + " = 24"

    # Split the extracted_text into left and right parts
    equation_part_list = extracted_text.split("=")
    if len(equation_part_list) == 2:
        left_side = equation_part_list[0]
        right_side = equation_part_list[1]
    elif len(equation_part_list) > 2:
        left_side = equation_part_list[0]
        right_side = equation_part_list[-1]
    #left_side, right_side = extracted_text.split("=")

    # Evaluate the right side
    try:
        right_value = eval(right_side.strip())
        if right_value != 24:
            return False
    except:
        return False

    # Prepare the number list as a multiset
    number_multiset = sorted(number_list)

    # Extract numbers and operators from the left side using regex
    # Supports LaTeX expressions as well
    left_side_numbers = re.findall(r'\d+', left_side)
    left_side_numbers = list(map(int, left_side_numbers))
    left_side_numbers_sorted = sorted(left_side_numbers)

    if number_multiset != left_side_numbers_sorted:
        return False

    # Replace LaTeX syntax with Python syntax for evaluation
    left_side = left_side.replace(r'\times', '*').replace(r'\div', '/')
    left_side = re.sub(r'\\frac{(\d+)}{(\d+)}', r'(\1/\2)', left_side)

    # Evaluate the left side
    try:
        left_value = eval(left_side.strip())
        if left_value != 24:
            return False
    except:
        return False
    return True

##### Letters #####
def evaluate_response_letters(word, target_letter, llm_response):
    """
    Evaluate the LLM's response for correctness.

    :param word: The test word
    :param target_letter: The letter that was counted
    :param llm_response: The response from the LLM
    :return: A tuple (is_correct, explanation)
    """
    # Extract count and positions from LLM response
    match = re.search(r"Count: (\d+), Positions: \[([\d, ]+)\]", llm_response)
    if not match:
        return False, "Response format is incorrect"

    llm_count = int(match.group(1))
    llm_positions = [int(pos) for pos in match.group(2).split(',')]

    # Calculate correct count and positions
    correct_count = word.count(target_letter)
    correct_positions = [i + 1 for i, letter in enumerate(word) if letter == target_letter]

    if llm_count != correct_count:
        return False, f"Incorrect count. Expected {correct_count}, got {llm_count}"

    if set(llm_positions) != set(correct_positions):
        return False, f"Incorrect positions. Expected {correct_positions}, got {llm_positions}"

    return True, "Correct response"

def extract_equation_with_GPT4_letters(response):
    prompt = 'Your task is to extract the final answer from the given answer by another LLM:\n' \
             'Note that the final answer should follow strictly the format like Count: 5, Positions: [2, 4, 13, 17, 22], Count: 1, Positions: [5], ' \
             'Count: 4, Positions: [3, 11, 18, 24] \n' \
             'Here is the response, return your answer with the format <<<final answer>>>, like <<<Count: 4, Positions: [3, 11, 18, 24]>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer are answered, then directly answer <<<Count: 0, Positions: []>>>.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs=False)
    return extract_equation

def create_prompt_letters(word, target_letter):
    prompt = f"How many '{target_letter}'s are in the word '{word}' and what are their positions? The position in the word counts from 1 not 0.\n"
    prompt += "Surround the answer with <<<content>>>. "
    prompt += f"Please respond in the format: <<<Count: X, Positions: [Y, Z, ...]>>>, such as <<<Count: 2, Positions: [1, 3]>>>.\n" \
              f"Your answer:\n"
    return prompt

def read_words_from_file_letters(filename):
    """
    Read a list of words from a JSON file.

    :param filename: Name of the file to read from
    :return: List of words
    """
    with open(filename, 'r') as f:
        words = json.load(f)
    #print(f"Words read from {filename}")
    return words

#####Number Multipy#######
def read_dataset_number_multipy(dataset_dir: str) -> List[Dict]:
    """Read all puzzle instances from the dataset directory in sample_i order"""
    puzzles = []
    sample_id = 0

    for digit_num_list in [[1, 2, 4], [1, 3, 4], [1, 2, 2, 4]]:
        dir_digit_name = f'digit'
        for digit_num in digit_num_list:
            dataset_base_dir = os.path.join(dataset_dir, f'{dir_digit_name}')
            for i in range(0, 20):
                dataset_base_dir_sample = os.path.join(dataset_base_dir, f'sample_{i}')
                generated_num_list = read_value_list(dataset_base_dir_sample + f"/input_value_list.txt")
                target_answer = read_answer(dataset_base_dir_sample + f"/target_answer.txt")

                equation_prompt = f'{generated_num_list[0]}'
                for generated_num in generated_num_list[1:]:
                    equation_prompt += f'*{generated_num}'
                question = f'What is the result of ' + equation_prompt + '?'
                puzzles.append({
                    'digit_num': digit_num,
                    'sample_id': sample_id,
                    'question': question,
                    'solution_data': target_answer
                })
    return puzzles

def extract_equation_with_GPT4_number_multiply(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<number>>>, like <<<43243.4>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<No answer found>>>.\n' \
             'If there is equation in the answer but no final numbers, do not calculate the number by yourself.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs=False)
    return extract_equation

def validate_solution_number_multiply(llm_response, target_answer):
    if int(llm_response) == target_answer:
        return True
    return False
def format_number_with_commas(number):
    return f"{number:,}"
def read_value_list(file_path):
    with open(file_path, 'r') as f:
        value_list = f.read()
    return ast.literal_eval(value_list)

def read_answer(file_path):
    with open(file_path, 'r') as f:
        answer = f.read()
    return int(answer)

#####Big Bench Hard#######
def read_dataset_big_bench_hard(dataset_dir: str) -> List[Dict]:
    puzzles = []
    for env_name in ['date_understanding', 'web_of_lies',
                   'logical_deduction_seven_objects', 'navigate']:
        DATA_PATH = dataset_dir + f'/{env_name}.json'
        question_json_list = []
        with open(DATA_PATH, 'r') as file:
            for line in file:
                question_json_list.append(json.loads(line))

        for i in range(0, len(question_json_list[0]['examples']), 4):
            print(f'Sample num: {i} in {env_name}, total is: {len(question_json_list[0]["examples"])}')

            data = question_json_list[0]['examples'][i]
            question = data['input'] + f'\n' + f'\nOutput final answer with the format <<<answer>>>.'
            target_answer = data['target']

            puzzles.append({
                # 'digit_num': digit_num,
                # 'sample_id': sample_id,
                'example_data': data,
                'question': question,
                'solution_data': target_answer
            })

    return puzzles
def is_equiv_func_big_bench_hard(question, target_answer, response):
    input_prompt_equiv_func = r'Evaluate whether the answer given by another LLM is correct. ' \
                              r'I will give you the question, the answer from another LLM, and the target answer. ' \
                              r'Then you need to extract the final answer from the tested LLM response and evaluate the correctness of the answer. ' \
                              r'In the end of your response, answer with the list of both extracted answer and Correct/Wrong judgement surrounded by <<<>>>. ' \
                              r'The examples are: <<<[No clear answer, Wrong]>>>, <<<[Yes, Wrong]>>>, <<<[No, Wrong]>>>, <<<[(C) 03/07/2017, Correct]>>>, <<<[(F) 01/24/1947, Wrong]>>>, ' \
                              r'<<<[(D) The motorcyle is the oldest, Wrong]>>>, <<<[(F) The quail is the fourth from the left, Correct]>>>, ' \
                              r'<<<[True, Correct]>>>, <<<[True, Wrong]>>>, <<<[False, Correct]>>>. ' \
                              f'\nMost of the time, the final answer is displayed at the end of the LLM response. However, if there is no clear answer, just answer <<<[No clear answer, Wrong]>>>. ' \
                              f'Do not change or summarize the answer by yourself. Just compare tested LLM answer with the target answer. Especially whether the options are matching! ' \
                              f'For instance, A and D options are not matching, which is judged as wrong! ' \
                              f'Now the question is: {question}; the tested answer is: {response}; the target answer is: {target_answer}. Your evaluation answer: '

    response = GPT_response('', input_prompt_equiv_func, model_name='gpt-4o',
                            code_interpreter=False, user_prompt_list = [input_prompt_equiv_func], response_total_list = [], logprobs = False)
    return response



#####Big Bench Hard#######
def read_dataset_gsm(dataset_dir: str) -> List[Dict]:
    puzzles = []
    question_json_list = []
    with open(dataset_dir, 'r') as file:
        for line in file:
            question_json_list.append(json.loads(line))
    for i in range(1, len(question_json_list), 25):
        data = question_json_list[i]
        target_answer = data['target']
        question = data[
                       'input'] + f'\n' + f'\nOutput final answer with the format <<<answer>>> such as <<<123.42>>>, <<<125.0>>>, <<<-9867>>>.\nYour answer: '
        puzzles.append({
            'example_data': data,
            'question': question,
            'solution_data': target_answer
        })

    return puzzles
def extract_equation_with_GPT4_gsm(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<list>>>, like <<<43243.4>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<No answer found>>>.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def is_equiv_func_gsm(target_answer, extracted_text):
    input_prompt_equiv_func = r'Evaluate whether the following numerical pair has the same values.' \
                              r'Neglect the format difference and the extra text like units and names and equations.' \
                              r'The value can be regarded as the same if they are < 1e-3 relative difference.' \
                              r'The examples are: ("12", "12.0", True), ("5*sqrt(13)", "15.97112779602377", False),' \
                              r'("10\text{ inches}", "10.0", True), ("42", "41.99999999999998", True), ("frac{63}{64}", "0.984375", True),' \
                              r'("frac{5\sqrt{5}}{3}", "5\sqrt{5}/3", True), (\tfrac34, "3/4", True), ("frac{1033}{4}+30\sqrt{3}", "169.0", False), ("AB=12+12\sqrt{3}", "12(\sqrt{3} + 1)", True),' \
                              r'((18, -18), (18, -18), True). ' \
                              r'In the end of your response, answer <<<True>>> or <<<False>>>'
    input_prompt_equiv_func = input_prompt_equiv_func + f'\n({target_answer}, {extracted_text}), Your answer:'
    response = GPT_response('Your are a helpful checker for math expressions.', input_prompt_equiv_func, model_name='gpt-4o',
                            code_interpreter=False, user_prompt_list = [input_prompt_equiv_func], response_total_list = [], logprobs = False)
    return response

#####Math Geometry#######
def read_dataset_math_geometry(dataset_dir: str) -> List[Dict]:
    puzzles = []
    for i in range(0, 1000):
        problem_path = dataset_dir + f'/{i}.json'
        if os.path.exists(problem_path):
            with open(problem_path, 'r') as file:
                data = json.load(file)
            question = data[
                           'problem'] + f'\n' + f'\nOutput final answer with the format <<<answer>>> such as <<<123.42>>>, <<<125.0>>>, <<<-9867>>>.\nYour answer: '
            target = data['solution']
            target_answer = remove_boxed(last_boxed_only_string(target))
            puzzles.append({
                'example_data': data,
                'question': question,
                'solution_data': target_answer
            })

    return puzzles

def last_boxed_only_string(string):
    idx = string.rfind("\\boxed")
    if idx < 0:
        idx = string.rfind("\\fbox")
        if idx < 0:
            return None

    i = idx
    right_brace_idx = None
    num_left_braces_open = 0
    while i < len(string):
        if string[i] == "{":
            num_left_braces_open += 1
        if string[i] == "}":
            num_left_braces_open -= 1
            if num_left_braces_open == 0:
                right_brace_idx = i
                break
        i += 1

    if right_brace_idx == None:
        retval = None
    else:
        retval = string[idx:right_brace_idx + 1]

    return retval

def remove_boxed(s):
    left = "\\boxed{"
    try:
        assert s[:len(left)] == left
        assert s[-1] == "}"
        return s[len(left):-1]
    except:
        return None

def extract_equation_with_GPT4_math_geometry(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<list>>>, like <<<43243.4>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<No answer found>>>.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def is_equiv_func_math_geometry(target_answer, extracted_text):
    input_prompt_equiv_func = r'Evaluate whether the following numerical pair has the same values.' \
                              r'Neglect the format difference and the extra text like units and names and equations.' \
                              r'The value can be regarded as the same if they are < 1e-3 relative difference.' \
                              r'The examples are: ("12", "12.0", True), ("5*sqrt(13)", "15.97112779602377", False),' \
                              r'("10\text{ inches}", "10.0", True), ("42", "41.99999999999998", True), ("frac{63}{64}", "0.984375", True),' \
                              r'("frac{5\sqrt{5}}{3}", "5\sqrt{5}/3", True), (\tfrac34, "3/4", True), ("frac{1033}{4}+30\sqrt{3}", "169.0", False), ("AB=12+12\sqrt{3}", "12(\sqrt{3} + 1)", True),' \
                              r'((18, -18), (18, -18), True). ' \
                              r'In the end of your response, answer <<<True>>> or <<<False>>>'
    input_prompt_equiv_func = input_prompt_equiv_func + f'\n({target_answer}, {extracted_text}), Your answer:'
    response = GPT_response('Your are a helpful checker for math expressions.', input_prompt_equiv_func, model_name='gpt-4o',
                            code_interpreter=False, user_prompt_list = [input_prompt_equiv_func], response_total_list = [], logprobs = False)
    return response

#####Math Counting and Probability#######
def read_dataset_math_counting_and_probability(dataset_dir: str) -> List[Dict]:
    puzzles = []
    for i in range(0, 1000):
        problem_path = dataset_dir + f'/{i}.json'
        if os.path.exists(problem_path):
            with open(problem_path, 'r') as file:
                data = json.load(file)
            question = data[
                           'problem'] + f'\n' + f'\nOutput final answer with the format <<<answer>>> such as <<<123.42>>>, <<<125.0>>>, <<<-9867>>>.\nYour answer: '
            target = data['solution']
            target_answer = remove_boxed(last_boxed_only_string(target))
            puzzles.append({
                'example_data': data,
                'question': question,
                'solution_data': target_answer
            })
    return puzzles

def extract_equation_with_GPT4_math_counting_and_probability(response):
    prompt = 'Your task is to extract the final numerical answer of the given answer by another LLM:\n' \
             'Here is the response, return your answer with the format <<<list>>>, like <<<43243.4>>>.\n' \
             'If the input text does not have <<<>>> and is already the pure answer, add <<<>>> and return your answer.\n' \
             'Note that if you find no final answer is answered, then directly answer <<<No answer found>>>.\n' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list = [prompt + response], response_total_list = [], logprobs = False)
    return extract_equation

def is_equiv_func_math_counting_and_probability(target_answer, extracted_text):
    input_prompt_equiv_func = r'Evaluate whether the following numerical pair has the same values.' \
                              r'Neglect the format difference and the extra text like units and names and equations.' \
                              r'The value can be regarded as the same if they are < 1e-3 relative difference.' \
                              r'The examples are: ("12", "12.0", True), ("5*sqrt(13)", "15.97112779602377", False),' \
                              r'("10\text{ inches}", "10.0", True), ("42", "41.99999999999998", True), ("frac{63}{64}", "0.984375", True),' \
                              r'("frac{5\sqrt{5}}{3}", "5\sqrt{5}/3", True), (\tfrac34, "3/4", True), ("frac{1033}{4}+30\sqrt{3}", "169.0", False), ("AB=12+12\sqrt{3}", "12(\sqrt{3} + 1)", True),' \
                              r'((18, -18), (18, -18), True). ' \
                              r'In the end of your response, answer <<<True>>> or <<<False>>>'
    input_prompt_equiv_func = input_prompt_equiv_func + f'\n({target_answer}, {extracted_text}), Your answer:'
    response = GPT_response('Your are a helpful checker for math expressions.', input_prompt_equiv_func, model_name='gpt-4o',
                            code_interpreter=False, user_prompt_list = [input_prompt_equiv_func], response_total_list = [], logprobs = False)
    return response

######BoxNet1######
def read_dataset_boxnet1(dataset_dir: str) -> List[Dict]:
    puzzles = []
    for pg_row_num, pg_column_num in [(1, 2), (2, 2), (2, 4)]:
        for iteration_num in range(10):
            print(f'Row num is: {pg_row_num}, Column num is: {pg_column_num}, Iteration num is: {iteration_num}\n\n')
            with open(
                    dataset_dir + f'/env_pg_state_{pg_row_num}_{pg_column_num}/pg_state{iteration_num}/pg_state{iteration_num}.json',
                    'r') as file:
                pg_dict = json.load(file)
            file.close()

            pg_dict_initial = copy.deepcopy(pg_dict)
            prompt = create_prompt(pg_row_num, pg_column_num, pg_dict)
            question = prompt
            puzzles.append({
                'pg_dict': pg_dict,
                'pg_dict_initial': pg_dict_initial,
                'question': question
            })
    return puzzles

def create_prompt(pg_row_num, pg_column_num, pg_dict):
    state_update_prompt = state_update_func(pg_row_num, pg_column_num, pg_dict)
    prompt = '''
You are a central planner tasked with directing agents in a grid-like field to move colored boxes to their corresponding color-coded targets. Each agent occupies a 1x1 square and can only interact with objects within its square. Agents can move a box to an adjacent square or directly to a target square of the same color. A square may contain multiple boxes and targets.



The squares are identified by their center coordinates (e.g., square[0.5, 0.5]). Actions are formatted as: move(box_color, destination), where box_color is the color of the box and destination is either a target of the same color or an adjacent square.



Your objective is to create a sequence of action plans that instructs each agent to match all boxes to their color-coded targets in the most efficient manner.



Please adhere to the following rules when specifying your action plan:



1. **Single Action per Agent**: Assign only one action to each agent at a time. However, the final answer shoule be a list of action plans for multiple steps.



2. **Unique Agent Keys**: Use unique keys for each agent in the JSON format action plan. The key should be the agent's coordinates in the format "Agent[x, y]".



3. **Prioritize Matching Boxes to Targets**: Always prioritize actions that will match a box to its target over moving a box to an adjacent square.



4. **Sequential Action Planning**: The whole returned answer should be a list of action plans for multiple steps, do not just return one step plan.



5. **Clear Formatting**: Ensure the action plan is clearly formatted in JSON, with each agent's action specified as a key-value pair.



6. **Conflict Resolution**: Ensure that no two agents are assigned actions that would interfere with each other.



7. **Optimize Efficiency**: Aim to minimize the number of moves required to match all boxes with their targets.



Here is the format for your action plan:
```json
[{"Agent[0.5, 0.5]":"move(box_blue, square[0.5, 1.5])", "Agent[1.5, 0.5]":"move(box_red, target_red)"}, {"Agent[0.5, 1.5]":"move(box_blue, target_blue)", "Agent[2.5, 0.5]":"move...}, {...}...]
```
Include an agent in the action plan only if it has a task to perform next.

'''
    prompt += "Surround the answer with <<<content>>>. \n"
    prompt += f'''
    The current left boxes and agents are:
    {state_update_prompt}\n
    '''
    prompt += '''
    Please respond in the format: <<<list of action dictionary>>>, such as <<<[{"Agent[0.5, 0.5]":"move(box_blue, square[0.5, 1.5])", "Agent[1.5, 0.5]":"move...}, {"Agent[0.5, 1.5]":"move(box_blue, target_blue)"}, {...}...]>>>.\n
    Your answer:\n
    '''
    return prompt

def score_in_training_set(pg_dict, response):
    success_failure = ''
    try:
        original_response_dict_list = json.loads(response)
        for original_response_dict in original_response_dict_list:
            for key, value in original_response_dict.items():
                coordinates = tuple(map(float, re.findall(r"\d+\.?\d*", key)))
                # match the item and location in the value
                match = re.match(r"move\((.*?),\s(.*?)\)", value)
    except:
        success_failure = 'response in the wrong format'

    if success_failure == 'response in the wrong format':
        print('\nResponse in the wrong format!\n')
        return pg_dict, success_failure
    elif success_failure == '':
        pg_dict_returned = action_from_response(pg_dict, original_response_dict_list)
        count = 0
        for key, value in pg_dict_returned.items():
            count += len(value)
        if count == 0:
            success_failure = 'success'
        elif success_failure == '':
            success_failure = 'failure after full execution'
        return pg_dict_returned, success_failure
def action_from_response(pg_dict_input, original_response_dict_list):
    pg_dict_current = copy.deepcopy(pg_dict_input)

    for original_response_dict in original_response_dict_list:
        transformed_dict = {}
        for key, value in original_response_dict.items():
            coordinates = tuple(map(float, re.findall(r"\d+\.?\d*", key)))
            match = re.match(r"move\((.*?),\s(.*?)\)", value)
            if match:
                item, location = match.groups()
                if "square" in location:
                    location = tuple(map(float, re.findall(r"\d+\.?\d*", location)))
                transformed_dict[coordinates] = [item, location]

        # Process each move with the current state
        for key, value in transformed_dict.items():
            current_pos = f"{key[0]}_{key[1]}"

            # Check if this is a box-target matching move
            if (value[0] in pg_dict_current[current_pos] and
                    isinstance(value[1], str) and
                    value[1] in pg_dict_current[current_pos] and
                    value[0].startswith('box_') and
                    value[1].startswith('target_') and
                    value[0][4:] == value[1][7:]):
                # Remove both box and target when matched
                pg_dict_current[current_pos].remove(value[0])
                pg_dict_current[current_pos].remove(value[1])

            # Check if this is a movement to another square
            elif (value[0] in pg_dict_current[current_pos] and
                  isinstance(value[1], tuple)):  # Only check coordinates for square movements
                # Calculate if move is to adjacent square
                if ((np.abs(key[0] - value[1][0]) == 0 and np.abs(key[1] - value[1][1]) == 1) or
                        (np.abs(key[0] - value[1][0]) == 1 and np.abs(key[1] - value[1][1]) == 0)):
                    # Move box to new location
                    target_pos = f"{value[1][0]}_{value[1][1]}"
                    pg_dict_current[current_pos].remove(value[0])
                    pg_dict_current[target_pos].append(value[0])

    return pg_dict_current

def extract_equation_with_GPT4_boxnet1(response):
    prompt = 'Your task is to extract the final answer from the given answer by another LLM:\n' \
             'Note that the equation should be in the form like <<<answer>>>, <<<[{"Agent[0.5, 0.5]":"move(box_blue, square[0.5, 1.5])", "Agent[1.5, 0.5]":"move...}, {"Agent[0.5, 1.5]":"move(box_blue, target_blue)}, {...}...]>>>, \n' \
             'Here is the reponse, return your answer with the format <<<equation>>>, like <<<[{"Agent[0.5, 0.5]":"move(box_blue, square[0.5, 1.5])", "Agent[1.5, 0.5]":"move...}, {"Agent[0.5, 1.5]":"move(box_blue, target_blue)}, {...}...]>>>. ' \
             'Input text: ' \

    extract_equation = GPT_response('', prompt + response, model_name='gpt-4o', code_interpreter=False, user_prompt_list=[prompt + response], response_total_list=[], logprobs=False)
    return extract_equation
