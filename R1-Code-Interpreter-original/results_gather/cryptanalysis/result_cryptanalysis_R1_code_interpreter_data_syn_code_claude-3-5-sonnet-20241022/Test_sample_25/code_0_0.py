def check_numbers(guess_nums, actual_nums, feedback):
    correct_pos = 0
    correct_val = 0
    too_small = 0
    
    for i in range(2):
        if guess_nums[i] == actual_nums[i]:
            correct_pos += 1
        elif guess_nums[i] in actual_nums:
            correct_val += 1
        elif guess_nums[i] < actual_nums[i]:
            too_small += 1
    
    return correct_pos, correct_val, too_small

def check_letters(guess_letters, actual_letters, feedback):
    correct_pos = 0
    correct_val = 0
    
    for i in range(2):
        if guess_letters[i] == actual_letters[i]:
            correct_pos += 1
        elif guess_letters[i] in actual_letters:
            correct_val += 1
    
    return correct_pos, correct_val

def matches_feedback(guess, actual, feedback_type):
    if feedback_type == 1:  # 91SF
        num_correct_pos, num_correct_val, _ = check_numbers(guess[:2], actual[:2], feedback_type)
        let_correct_pos, let_correct_val = check_letters(guess[2:], actual[2:], feedback_type)
        return (num_correct_pos == 0 and 
                let_correct_pos == 1 and let_correct_val == 0)
    
    elif feedback_type == 2:  # 03HF
        num_correct_pos, num_correct_val, num_too_small = check_numbers(guess[:2], actual[:2], feedback_type)
        let_correct_pos, let_correct_val = check_letters(guess[2:], actual[2:], feedback_type)
        return (num_correct_pos == 1 and num_too_small == 1 and 
                let_correct_pos == 0 and let_correct_val == 0)
    
    elif feedback_type == 3:  # 75CP
        num_correct_pos, num_correct_val, _ = check_numbers(guess[:2], actual[:2], feedback_type)
        let_correct_pos, let_correct_val = check_letters(guess[2:], actual[2:], feedback_type)
        return (num_correct_pos == 0 and 
                let_correct_pos == 0 and let_correct_val == 1)
    
    elif feedback_type == 4:  # 64YB
        num_correct_pos, num_correct_val, _ = check_numbers(guess[:2], actual[:2], feedback_type)
        let_correct_pos, let_correct_val = check_letters(guess[2:], actual[2:], feedback_type)
        return (num_correct_pos == 1 and num_correct_val == 0 and 
                let_correct_pos == 0 and let_correct_val == 0)

import string

# Generate all possible combinations
numbers = list(range(10))
letters = list(string.ascii_uppercase)

guesses = [
    ('91SF', 1),
    ('03HF', 2),
    ('75CP', 3),
    ('64YB', 4)
]

solutions = []

for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:
                    continue
                    
                candidate = f"{n1}{n2}{l1}{l2}"
                valid = True
                
                for guess, feedback_type in guesses:
                    if not matches_feedback(guess, candidate, feedback_type):
                        valid = False
                        break
                
                if valid:
                    solutions.append([str(n1), str(n2), l1, l2])

print(solutions)