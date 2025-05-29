from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = ['0', '3', '6', '8']
possible_letters = ['Q', 'X', 'Z']

# Feedback rules
feedback = [
    ('83AJ', (1, 1, 0, 0)),  # (correct number, wrong position, incorrect number too large, incorrect letter too early)
    ('72XZ', (0, 0, 1, 1)),  # (correct letter, correct position, incorrect letter, incorrect)
    ('03LZ', (1, 1, 0, 1)),  # (correct number, wrong position, correct letter, correct position)
    ('25KR', (0, 0, 0, 0)),  # (all incorrect)
    ('36TF', (1, 1, 0, 0)),  # (correct number, correct position, incorrect number too large, incorrect letters)
    ('15JN', (0, 0, 0, 0)),  # (all incorrect)
    ('16FQ', (0, 0, 1, 1))   # (correct letter, wrong position, incorrect letter too early)
]

def is_valid_combination(combo):
    for guess, (cn, wp, il, ie) in feedback:
        numbers, letters = combo[:2], combo[2:]
        guess_numbers, guess_letters = guess[:2], guess[2:]
        
        # Check numbers
        correct_numbers = sum(n in numbers for n in guess_numbers)
        wrong_position_numbers = sum(n in numbers and numbers.index(n) != guess_numbers.index(n) for n in guess_numbers)
        incorrect_large_numbers = sum(n > max(numbers) for n in guess_numbers)
        
        # Check letters
        correct_letters = sum(l in letters for l in guess_letters)
        correct_position_letters = sum(l in letters and letters.index(l) == guess_letters.index(l) for l in guess_letters)
        incorrect_early_letters = sum(l < min(letters) for l in guess_letters)
        
        if not (correct_numbers == cn and wrong_position_numbers == wp and
                incorrect_large_numbers == il and correct_letters == correct_position_letters and
                incorrect_early_letters == ie):
            return False
    return True

# Generate all permutations of possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for let_perm in permutations(possible_letters, 2):
        combination = num_perm + let_perm
        if is_valid_combination(combination):
            print(f"<<< {list(combination)} >>>")
            break