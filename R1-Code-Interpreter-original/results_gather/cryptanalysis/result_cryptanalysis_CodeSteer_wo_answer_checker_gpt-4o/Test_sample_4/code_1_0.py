from itertools import permutations

# Define possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Constraints based on feedback
def is_valid_combination(num1, num2, let1, let2):
    # Guess 1: 41CK
    if not ((num1 == '4' and num2 != '1') or (num2 == '4' and num1 != '1')):
        return False
    if not (let1 == 'C' or let2 == 'C'):
        return False
    if let1 == 'K' or let2 == 'K':
        return False

    # Guess 2: 64DE
    if not ((num1 == '6' and num2 != '4') or (num2 == '6' and num1 != '4')):
        return False
    if let1 in 'DE' or let2 in 'DE':
        return False

    # Guess 3: 87JY
    if num1 in '87' or num2 in '87':
        return False
    if let1 in 'JY' or let2 in 'JY':
        return False

    # Guess 4: 12OD
    if num1 in '12' or num2 in '12':
        return False
    if not (let1 == 'D' or let2 == 'D'):
        return False
    if let1 == 'O' or let2 == 'O':
        return False

    return True

# Generate all possible combinations
for num1, num2 in permutations(numbers, 2):
    for let1, let2 in permutations(letters, 2):
        if is_valid_combination(num1, num2, let1, let2):
            password = [num1, num2, let1, let2]
            print(f"<<< {password} >>>")
            break