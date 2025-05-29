from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = ['9', '0']
possible_letters = ['X', 'Y', 'Z']

# All possible combinations of two numbers and two letters
combinations = permutations(possible_numbers, 2)
letter_combinations = permutations(possible_letters, 2)

# Function to check if a combination satisfies all conditions
def is_valid_combination(numbers, letters):
    # Guess 18WV
    if not ((numbers[0] == '1' or numbers[1] == '8') and (letters[0] == 'W' or letters[1] == 'V')):
        return False
    # Guess 76GH
    if not (numbers[0] > '7' and numbers[1] > '6' and letters[0] > 'H' and letters[1] > 'H'):
        return False
    # Guess 75LE
    if not (numbers[0] > '7' and numbers[1] > '5' and letters[0] > 'L' and letters[1] > 'E'):
        return False
    # Guess 12BZ
    if not ((numbers[0] > '2' and numbers[1] > '1') and (letters[0] == 'B' or letters[1] == 'Z')):
        return False
    return True

# Iterate over all combinations
for num_comb in combinations:
    for let_comb in letter_combinations:
        if is_valid_combination(num_comb, let_comb):
            password = list(num_comb + let_comb)
            print(f"<<< {password} >>>")
            break