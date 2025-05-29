from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
conditions = [
    ('80KR', (1, 1, 0, 0)),  # (correct number wrong position, incorrect number too small, incorrect letters)
    ('87NP', (1, 1, 1, 0)),  # (correct number wrong position, incorrect number, correct letter correct position, incorrect letter)
    ('47AZ', (0, 0, 1, 1)),  # (incorrect numbers, correct letter correct position, incorrect letter too late)
    ('15RF', (0, 0, 0, 0))   # (incorrect numbers too small, incorrect letters)
]

# Function to check a combination against the conditions
def check_combination(combo):
    for guess, feedback in conditions:
        num_correct_wrong_pos = sum((c in guess[:2]) and (c != guess[i]) for i, c in enumerate(combo[:2]))
        num_incorrect_too_small = sum((c not in guess[:2]) and (c < guess[i]) for i, c in enumerate(combo[:2]))
        letter_correct_pos = sum(c == guess[i] for i, c in enumerate(combo[2:]))
        letter_incorrect = sum(c not in guess[2:] for c in combo[2:])
        
        if (num_correct_wrong_pos != feedback[0] or
            num_incorrect_too_small != feedback[1] or
            letter_correct_pos != feedback[2] or
            letter_incorrect != feedback[3]):
            return False
    return True

# Generate all possible combinations of two numbers and two letters
for num_combo in permutations(numbers, 2):
    for letter_combo in permutations(letters, 2):
        combo = num_combo + letter_combo
        if check_combination(combo):
            print(f"<<< {list(combo)} >>>")
            break