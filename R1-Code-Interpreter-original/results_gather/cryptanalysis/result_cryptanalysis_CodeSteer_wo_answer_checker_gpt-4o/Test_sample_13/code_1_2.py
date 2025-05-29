from itertools import permutations

# Define possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
feedback = [
    ('95VB', (1, 1, 0, 0)),  # (one number correct but wrong position, one number incorrect, both letters incorrect)
    ('34OH', (0, 2, 1, 1)),  # (both numbers incorrect, one letter correct and correct position, one letter incorrect and too early)
    ('63TK', (0, 2, 0, 2)),  # (both numbers incorrect, both letters incorrect)
    ('61XN', (1, 1, 1, 1))   # (one number correct but wrong position, one number incorrect, one letter correct and correct position, one letter incorrect and too late)
]

# Function to check if a combination satisfies the feedback
def satisfies_feedback(combo):
    for guess, (num_correct_wrong_pos, num_incorrect, letter_correct_pos, letter_incorrect) in feedback:
        # Check numbers
        num_correct_wrong_pos_actual = sum((c in guess[:2]) for c in combo[:2]) - sum((c == guess[i]) for i, c in enumerate(combo[:2]))
        num_incorrect_actual = 2 - sum((c in guess[:2]) for c in combo[:2])
        
        # Check letters
        letter_correct_pos_actual = sum((c == guess[i+2]) for i, c in enumerate(combo[2:]))
        letter_incorrect_actual = 2 - sum((c in guess[2:]) for c in combo[2:])
        
        # Check if the actual feedback matches the expected feedback
        if (num_correct_wrong_pos_actual != num_correct_wrong_pos or
            num_incorrect_actual != num_incorrect or
            letter_correct_pos_actual != letter_correct_pos or
            letter_incorrect_actual != letter_incorrect):
            return False
    return True

# Iterate over all possible combinations
for num_combo in permutations(numbers, 2):
    for letter_combo in permutations(letters, 2):
        combo = num_combo + letter_combo
        if satisfies_feedback(combo):
            print(f"<<< {list(combo)} >>>")