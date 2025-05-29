from itertools import permutations

# Possible numbers and letters based on constraints
possible_numbers = [str(i) for i in range(4, 7)]  # Numbers between 4 and 6
possible_letters = [chr(i) for i in range(ord('E'), ord('O') + 1)]  # Letters between E and O

# Function to check if a combination satisfies all conditions
def is_valid_combination(numbers, letters):
    # Guess 1: 03OQ
    if not (numbers[0] > '0' and numbers[1] > '3' and letters[0] < 'O' and letters[1] < 'Q'):
        return False
    # Guess 2: 76EC
    if not (numbers[0] < '7' and numbers[1] < '6' and letters[0] > 'E' and letters[1] > 'C'):
        return False
    # Guess 3: 83FO
    if not ((letters[0] == 'O' and letters[1] != 'F') or (letters[1] == 'O' and letters[0] != 'F')):
        return False
    # Guess 4: 47OF
    if not ((numbers[0] == '4' and numbers[1] < '7') and ((letters[0] == 'F' and letters[1] != 'O') or (letters[1] == 'F' and letters[0] != 'O'))):
        return False
    # Guess 5: 87LK
    if not ((letters[0] == 'K' and letters[1] != 'L') or (letters[1] == 'K' and letters[0] != 'L')):
        return False
    return True

# Find the correct combination
for num_comb in permutations(possible_numbers, 2):
    for let_comb in permutations(possible_letters, 2):
        if is_valid_combination(num_comb, let_comb):
            password = list(num_comb + let_comb)
            print(password)
            break