from itertools import permutations

def is_valid_combination(numbers, letters):
    # Check against Guess 95WC
    if not ((numbers[0] in '95' and numbers[0] != '9') or (numbers[1] in '95' and numbers[1] != '5')):
        return False
    if '9' in numbers and '5' in numbers:
        return False
    if letters[0] != 'W':
        return False
    if letters[1] <= 'C':
        return False

    # Check against Guess 12XN
    if not ((numbers[0] == '1' and numbers[1] != '2') or (numbers[1] == '2' and numbers[0] != '1')):
        return False
    if 'X' in letters or 'N' in letters:
        return False

    # Check against Guess 92LS
    if not ((numbers[0] in '92' and numbers[0] != '9') or (numbers[1] in '92' and numbers[1] != '2')):
        return False
    if '9' in numbers and '2' in numbers:
        return False
    if not ((letters[0] == 'L' and letters[1] != 'S') or (letters[1] == 'L' and letters[0] != 'S')):
        return False

    return True

# Possible numbers and letters
possible_numbers = '0123456789'
possible_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Generate all permutations of two numbers and two letters
for numbers in permutations(possible_numbers, 2):
    for letters in permutations(possible_letters, 2):
        if is_valid_combination(numbers, letters):
            password = list(numbers + letters)
            print(f"<<< {password} >>>")
            break