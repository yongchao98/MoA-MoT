from itertools import permutations

# Define possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
feedback = [
    ('95VB', (1, 0, 0)),  # (correct but wrong position, incorrect, incorrect)
    ('34OH', (0, 1, 1)),  # (incorrect, correct and correct position, incorrect)
    ('63TK', (0, 0, 0)),  # (incorrect, incorrect)
    ('61XN', (1, 0, 1))   # (correct but wrong position, incorrect, correct and correct position, incorrect)
]

# Function to check if a combination satisfies the feedback
def satisfies_feedback(combo):
    for guess, (num_feedback, letter_feedback) in feedback:
        num_correct_wrong_pos = sum((c in guess[:2]) for c in combo[:2]) - sum((c == guess[i]) for i, c in enumerate(combo[:2]))
        num_incorrect = 2 - sum((c in guess[:2]) for c in combo[:2])
        letter_correct_pos = sum((c == guess[i+2]) for i, c in enumerate(combo[2:]))
        letter_incorrect = 2 - sum((c in guess[2:]) for c in combo[2:])
        
        if (num_correct_wrong_pos != num_feedback[0] or
            num_incorrect != num_feedback[1] or
            letter_correct_pos != letter_feedback[0] or
            letter_incorrect != letter_feedback[1]):
            return False
    return True

# Iterate over all possible combinations
for num_combo in permutations(numbers, 2):
    for letter_combo in permutations(letters, 2):
        combo = num_combo + letter_combo
        if satisfies_feedback(combo):
            print(f"<<< {list(combo)} >>>")