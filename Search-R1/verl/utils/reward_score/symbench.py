import re
import string
import random
from r1_code_inter.Logic_Game_func import verify_solution_func_gather, load_task_dataset
from r1_code_inter.prompt import R1_code_interpreter_data_syn_prompt1

def compute_score_symbench(solution_str, ground_truth, method='strict', format_score=0., score=1.):
    """The scoring function for exact match (EM).

    Args:
        solution_str: the solution text
        ground_truth: the ground truth
        method: the method to extract the solution, choices are 'strict' and 'flexible'
        format_score: the score for the format
        score: the score for the correct answer
    """
    i = int(ground_truth['index']); task_name = ground_truth['task_name']; response = solution_str

    solution_list, question_list, target_list, puzzles, solution_data_list, question_constrained_list, question_matrix_list, number_list, word_list, letter_list, save_input_dir = load_task_dataset(
        task_name, '')

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

    if isinstance(response, (bytes, bytearray)):
        response = response.decode('utf-8', errors='ignore')

    if not isinstance(response, str):
        print(f"Error: response is not a string, but {type(response)}")
        return 0

    if type(response) != str:
        print(f"Error: response is not a string, but {type(response)}")
        return 0
    try:
        True_false_result_1, True_false_result_2 = verify_solution_func_gather(i, task_name, response, question,
                                                                           solution, target, puzzles, solution_data_list,
                                                                           solution_list, question_constrained_list, question_matrix_list, number_list_item, word, letter)
    except:
        print(f"Error: verify_solution_func_gather failed for task {task_name} with index {i}")
        return 0

    ### Remain unchanged
    if True_false_result_1 == False and True_false_result_2 == False:
        return 0
    else:
        return score

    '''
    answer = extract_solution(solution_str=solution_str)
    do_print = random.randint(1, 64) == 1
    if do_print:
        print(f"--------------------------------")
        print(f"Golden answers: {ground_truth['target']}")
        print(f"Extracted answer: {answer}")
        print(f"Solution string: {solution_str}")
    '''
