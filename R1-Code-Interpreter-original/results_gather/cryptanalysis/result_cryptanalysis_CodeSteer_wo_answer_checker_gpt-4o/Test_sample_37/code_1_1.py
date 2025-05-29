from itertools import permutations

# Possible numbers and letters based on constraints
possible_numbers = [str(i) for i in range(4, 7)]  # Numbers between 4 and 6
possible_letters = [chr(i) for i in range(ord('E'), ord('O') + 1)]  # Letters between E and O

# All possible combinations of two numbers and two letters
all_combinations = permutations(possible_numbers, 2) + permutations(possible_letters, 2)

# Function to check if a combination satisfies all conditions
def is_valid_combination(numbers, letters):
    # Check against each guess
    # Guess 1: 03OQ
    if not (numbers[0] > '0' and numbers[1] > '3' and letters[0] < 'O' and letters[1] < 'Q'):
        return False
    # Guess 2: 76EC
    if not (numbers[0] < '7' and numbers[1] < '6' and letters[0] > 'E' and letters[1] > 'C'):
        return False
    # Guess 3: 83FO
    if not (letters[0] == 'O' or letters[1] == 'O'):
        return False
    # Guess 4: 47OF
    if not (numbers[0] == '4' or numbers[1] < '7' and letters[0] == 'F' or letters[1] == 'F'):
        return False
    # Guess 5: 87LK
    if not (letters[0] == 'K' or letters[1] == 'K'):
        return False
    return True

# Find the correct combination
for num_comb in permutations(possible_numbers, 2):
    for let_comb in permutations(possible_letters, 2):
        if is_valid_combination(num_comb, let_comb):
            password = list(num_comb + let_comb)
            print(password)
            break