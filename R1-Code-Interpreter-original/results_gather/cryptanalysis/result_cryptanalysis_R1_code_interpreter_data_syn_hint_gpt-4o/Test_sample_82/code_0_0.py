from itertools import permutations

# Possible numbers and letters based on deductions
possible_numbers = ['0', '1', '6', '9']
possible_letters = ['E', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R']

# Function to check if a guess matches the feedback
def check_guess(guess, feedback):
    num_correct = 0
    num_wrong_position = 0
    for i, (g, f) in enumerate(zip(guess, feedback)):
        if g == f:
            num_correct += 1
        elif g in feedback:
            num_wrong_position += 1
    return num_correct, num_wrong_position

# Iterate over all permutations of possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        combination = num_perm + letter_perm
        # Check against all feedback
        if (check_guess(combination, '87CF') == (0, 0) and
            check_guess(combination, '71SM') == (0, 2) and
            check_guess(combination, '23AY') == (0, 0) and
            check_guess(combination, '53PD') == (0, 1) and
            check_guess(combination, '01UF') == (0, 1) and
            check_guess(combination, '43US') == (0, 0)):
            print(combination)
            break