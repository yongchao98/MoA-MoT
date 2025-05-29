from itertools import permutations

# Define possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
conditions = [
    ('08VY', (1, 1, 0, 0)),  # Guess 1
    ('07CX', (1, 1, 0, 0)),  # Guess 2
    ('45NJ', (0, 0, 1, 0)),  # Guess 3
    ('28YQ', (0, 1, 0, 1))   # Guess 4
]

def check_combination(combo):
    for guess, feedback in conditions:
        num_correct_pos = sum(combo[i] == guess[i] for i in range(4))
        num_correct = sum(c in guess for c in combo) - num_correct_pos
        if (num_correct_pos, num_correct) != feedback:
            return False
    return True

# Generate all permutations of 2 numbers and 2 letters
for num_perm in permutations(numbers, 2):
    for letter_perm in permutations(letters, 2):
        combination = num_perm + letter_perm
        if check_combination(combination):
            print(f"<<< {list(combination)} >>>")
            break