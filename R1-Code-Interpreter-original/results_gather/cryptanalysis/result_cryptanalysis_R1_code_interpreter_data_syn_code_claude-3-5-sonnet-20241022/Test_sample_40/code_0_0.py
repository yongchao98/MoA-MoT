from itertools import permutations

def check_guess(candidate, guess, feedback):
    # Split feedback into numbers and letters
    num_correct_pos = 0
    num_wrong_pos = 0
    let_correct_pos = 0
    let_wrong_pos = 0
    
    # Check numbers (first two positions)
    for i in range(2):
        if candidate[i] == guess[i]:
            num_correct_pos += 1
        elif guess[i] in candidate[:2]:
            num_wrong_pos += 1
            
    # Check letters (last two positions)
    for i in range(2, 4):
        if candidate[i] == guess[i]:
            let_correct_pos += 1
        elif guess[i] in candidate[2:]:
            let_wrong_pos += 1
    
    return (num_correct_pos, num_wrong_pos, let_correct_pos, let_wrong_pos)

def matches_feedback(candidate, guess, expected):
    result = check_guess(candidate, guess, expected)
    return result == expected

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Store all guesses and their feedback
guesses = [
    ('17LN', (0, 0, 0, 0)),  # both numbers incorrect, both letters incorrect
    ('52MG', (1, 0, 0, 0)),  # one number correct position, one incorrect, both letters incorrect
    ('39MH', (0, 0, 0, 0)),  # both numbers incorrect, both letters incorrect
    ('85QP', (1, 0, 1, 0)),  # one number correct position, one incorrect, one letter correct position, one incorrect
    ('75OJ', (0, 0, 0, 0)),  # both numbers incorrect, both letters incorrect
    ('40AG', (0, 0, 0, 1)),  # both numbers incorrect, one letter wrong position, one incorrect
]

# Try all possible combinations
valid_combinations = []
for num_perm in permutations(numbers, 2):
    for let_perm in permutations(letters, 2):
        candidate = list(num_perm + let_perm)
        valid = True
        for guess, feedback in guesses:
            if not matches_feedback(candidate, list(guess), feedback):
                valid = False
                break
        if valid:
            valid_combinations.append(candidate)

print(valid_combinations)