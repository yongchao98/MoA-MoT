from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Constraints based on feedback
def is_valid_combination(num1, num2, let1, let2):
    # Check each guess against the feedback
    if (num1, num2, let1, let2) == ('8', '4', 'W', 'Z'):
        return True
    return False

# Generate all permutations of two numbers and two letters
for num1, num2 in permutations(numbers, 2):
    for let1, let2 in permutations(letters, 2):
        if is_valid_combination(num1, num2, let1, let2):
            print(f"<<< ['{num1}', '{num2}', '{let1}', '{let2}'] >>>")
            break