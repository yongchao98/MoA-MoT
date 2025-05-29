from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback from each guess
feedback = [
    ('17LN', (0, 0)),  # both numbers and letters are incorrect
    ('52MG', (1, 0)),  # one number correct and in the correct position, both letters incorrect
    ('39MH', (0, 0)),  # both numbers and letters are incorrect
    ('85QP', (1, 1)),  # one number correct and in the correct position, one letter correct and in the correct position
    ('75OJ', (0, 0)),  # both numbers and letters are incorrect
    ('40AG', (0, 1))   # both numbers incorrect, one letter correct but in the wrong position
]

# Function to check a combination against feedback
def check_combination(combo, feedback):
    for guess, (num_correct, let_correct) in feedback:
        num_match = sum(c1 == c2 for c1, c2 in zip(combo[:2], guess[:2]))
        let_match = sum(c1 == c2 for c1, c2 in zip(combo[2:], guess[2:]))
        if num_match != num_correct or let_match != let_correct:
            return False
    return True

# Generate all possible combinations of two numbers and two letters
for num_combo in permutations(numbers, 2):
    for let_combo in permutations(letters, 2):
        combo = ''.join(num_combo + let_combo)
        if check_combination(combo, feedback):
            password = list(combo)
            print(f"<<< {password} >>>")
            break