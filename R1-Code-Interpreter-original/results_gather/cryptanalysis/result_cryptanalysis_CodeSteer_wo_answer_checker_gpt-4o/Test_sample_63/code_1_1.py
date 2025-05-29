from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = [3, 6]
possible_letters = ['O', 'P']

# All possible combinations of two numbers and two letters
number_combinations = list(permutations(possible_numbers, 2))
letter_combinations = list(permutations(possible_letters, 2))

# Function to check if a combination satisfies all constraints
def is_valid_combination(numbers, letters):
    # Check against each guess and feedback
    # Guess 1: 38UW
    if not ((numbers[0] == 3 or numbers[1] == 3) and (numbers[0] != 8 and numbers[1] != 8)):
        return False
    # Guess 2: 69OG
    if not (numbers[0] == 6 and letters[0] == 'O'):
        return False
    # Guess 3: 83OP
    if not (numbers[1] == 3 and letters[0] == 'O'):
        return False
    # Guess 4: 48SP
    if not (numbers[0] != 4 and numbers[1] != 8):
        return False
    # Guess 5: 42ZE
    if not (numbers[0] != 4 and numbers[1] != 2):
        return False
    # Guess 6: 01WG
    if not (numbers[0] != 0 and numbers[1] != 1):
        return False
    # Guess 7: 07JX
    if not (numbers[0] != 0 and numbers[1] != 7):
        return False
    # Guess 8: 31HK
    if not ((numbers[0] == 3 or numbers[1] == 3) and (numbers[0] != 1 and numbers[1] != 1)):
        return False
    # Guess 9: 37VY
    if not ((numbers[0] == 3 or numbers[1] == 3) and (numbers[0] != 7 and numbers[1] != 7)):
        return False
    # Guess 10: 62QH
    if not (numbers[0] == 6 and numbers[1] != 2):
        return False
    # Guess 11: 94CL
    if not (numbers[0] != 9 and numbers[1] != 4):
        return False
    # Guess 12: 71TE
    if not (numbers[0] != 7 and numbers[1] != 1):
        return False
    return True

# Find the correct combination
for numbers in number_combinations:
    for letters in letter_combinations:
        if is_valid_combination(numbers, letters):
            password = [str(numbers[0]), str(numbers[1]), letters[0], letters[1]]
            print(f"<<< {password} >>>")
            break