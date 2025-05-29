from itertools import permutations

# Possible numbers and letters
numbers = ['0', '3']
letters = ['A', 'C', 'D', 'F', 'G', 'H', 'I', 'K', 'M', 'N', 'O', 'Q', 'S', 'T']

# Function to check if a combination is valid based on the conditions
def is_valid_combination(num1, num2, letter1, letter2):
    # Check against each condition
    # Condition 1: 85UJ
    if num1 in '85' or num2 in '85' or letter1 in 'UJ' or letter2 in 'UJ':
        return False
    # Condition 2: 19BP
    if num1 in '19' or num2 in '19':
        return False
    if (letter1 == 'B' and letter2 != 'P') or (letter2 == 'B' and letter1 != 'P'):
        return False
    # Condition 3: 78LU
    if (num1 == '7' and num2 != '8') or (num2 == '7' and num1 != '8'):
        return False
    if letter1 == 'L' and letter2 != 'U':
        return False
    # Condition 4: 61VZ
    if num1 in '61' or num2 in '61' or letter1 in 'VZ' or letter2 in 'VZ':
        return False
    # Condition 5: 41EL
    if num1 in '41' or num2 in '41':
        return False
    if (letter1 == 'E' and letter2 != 'L') or (letter2 == 'E' and letter1 != 'L'):
        return False
    # Condition 6: 12PR
    if num1 in '12' or num2 in '12' or letter1 in 'PR' or letter2 in 'PR':
        return False
    return True

# Generate all permutations of the numbers and letters
for num1, num2 in permutations(numbers, 2):
    for letter1, letter2 in permutations(letters, 2):
        if is_valid_combination(num1, num2, letter1, letter2):
            print([num1, num2, letter1, letter2])