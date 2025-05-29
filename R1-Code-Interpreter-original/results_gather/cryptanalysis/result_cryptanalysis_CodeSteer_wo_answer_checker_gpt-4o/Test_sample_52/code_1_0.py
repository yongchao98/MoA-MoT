from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Generate all possible combinations of two numbers and two letters
possible_combinations = [''.join(p) for p in permutations(numbers, 2)]
possible_combinations = [p + ''.join(l) for p in possible_combinations for l in permutations(letters, 2)]

# Feedback from each guess
feedback = [
    ('15NT', 'both numbers are incorrect; both letters are incorrect and too late in the alphabet.'),
    ('21OD', 'both numbers are incorrect; both letters are incorrect.'),
    ('51BI', 'both numbers are incorrect; both letters are incorrect.'),
    ('06GE', 'one number is correct and in the correct position; one number is incorrect and incorrect; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet.'),
    ('32FY', 'both numbers are incorrect; both letters are incorrect.'),
    ('13UN', 'both numbers are incorrect; both letters are incorrect and too late in the alphabet.'),
    ('37DU', 'both numbers are incorrect; both letters are incorrect.'),
    ('15IG', 'both numbers are incorrect; one letter is correct but in the wrong position; one letter is incorrect and too late in the alphabet.'),
    ('25IA', 'both numbers are incorrect; both letters are incorrect.'),
    ('93PA', 'one number is correct but in the wrong position; one number is incorrect and incorrect; both letters are incorrect.')
]

# Function to check if a combination is valid based on feedback
def is_valid_combination(combination, guess, feedback):
    num_correct = 0
    num_wrong_position = 0
    num_incorrect = 0

    # Check numbers
    for i in range(2):
        if combination[i] == guess[i]:
            num_correct += 1
        elif combination[i] in guess[:2]:
            num_wrong_position += 1
        else:
            num_incorrect += 1

    # Check letters
    for i in range(2, 4):
        if combination[i] == guess[i]:
            num_correct += 1
        elif combination[i] in guess[2:]:
            num_wrong_position += 1
        else:
            num_incorrect += 1

    # Validate based on feedback
    if 'both numbers are incorrect' in feedback and num_correct > 0:
        return False
    if 'one number is correct and in the correct position' in feedback and num_correct != 1:
        return False
    if 'one number is correct but in the wrong position' in feedback and num_wrong_position != 1:
        return False
    if 'both letters are incorrect' in feedback and num_correct > 0:
        return False
    if 'one letter is correct and in the correct position' in feedback and num_correct != 1:
        return False
    if 'one letter is correct but in the wrong position' in feedback and num_wrong_position != 1:
        return False

    return True

# Filter combinations based on feedback
for guess, fb in feedback:
    possible_combinations = [comb for comb in possible_combinations if is_valid_combination(comb, guess, fb)]

# Output the final valid combination
print(f"<<< {list(possible_combinations[0])} >>>")