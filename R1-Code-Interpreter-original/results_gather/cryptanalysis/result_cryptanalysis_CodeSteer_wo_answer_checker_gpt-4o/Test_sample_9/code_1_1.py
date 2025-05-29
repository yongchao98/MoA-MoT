from itertools import permutations

# Define the constraints based on the feedback
def is_valid_combination(numbers, letters):
    # Check against each guess
    # Guess 1: 64RY
    if set(numbers).intersection({'6', '4'}) or set(letters).intersection({'R', 'Y'}):
        return False
    # Guess 2: 51MB
    if set(numbers).intersection({'5', '1'}) or set(letters).intersection({'M', 'B'}):
        return False
    # Guess 3: 29PV
    if set(numbers).intersection({'2', '9'}) or set(letters).intersection({'P', 'V'}):
        return False
    # Guess 4: 71EW
    if '7' in numbers and '1' in numbers:
        return False
    if 'E' in letters and 'W' in letters:
        return False
    if '7' not in numbers and '1' not in numbers:
        return False
    if 'E' not in letters and 'W' not in letters:
        return False
    return True

# Possible numbers and letters
possible_numbers = {'0', '3', '7'}
possible_letters = {'E', 'W'}

# Generate all permutations of the possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for let_perm in permutations(possible_letters, 2):
        if is_valid_combination(num_perm, let_perm):
            password = list(num_perm) + list(let_perm)
            print(f"<<< {password} >>>")
            break