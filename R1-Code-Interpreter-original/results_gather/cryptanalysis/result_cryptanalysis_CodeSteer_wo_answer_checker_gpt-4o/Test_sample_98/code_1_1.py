from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Constraints based on feedback
def is_valid_combination(num1, num2, let1, let2):
    # Convert to strings for easier comparison
    num1, num2 = str(num1), str(num2)
    
    # Check each condition
    conditions = [
        (num1 != '7' and num2 != '5' and let1 > 'Q' and let2 > 'G'),  # Condition 1
        (num1 > '3' and num2 > '0' and let1 > 'K' and let2 > 'J'),    # Condition 2
        (num1 != '9' and num2 != '7' and let1 > 'K' and let2 > 'P'),  # Condition 3
        (num1 != '1' and num2 != '6' and let1 > 'K' and let2 > 'A'),  # Condition 4
        ((num1 == '8' or num2 == '8') and (num1 != '6' and num2 != '6') and (let1 == 'Y' or let2 == 'Y') and (let1 > 'C' and let2 > 'C')),  # Condition 5
        (num1 != '9' and num2 != '1' and let1 > 'M' and let2 > 'S'),  # Condition 6
        (num1 > '3' and num2 > '1' and let1 != 'W' and let2 != 'R'),  # Condition 7
        ((num1 == '8' or num2 == '8') and (num1 != '7' and num2 != '7') and let1 != 'Z' and let2 != 'L'),  # Condition 8
        ((num1 == '4' or num2 == '4') and (num1 != '5' and num2 != '5') and let1 != 'W' and let2 != 'B'),  # Condition 9
        (num1 == '4' and num2 == '8' and let1 != 'Q' and let2 != 'W'),  # Condition 10
        ((num1 == '8' or num2 == '8') and (num1 != '6' and num2 != '6') and let1 > 'U' and let2 > 'T')   # Condition 11
    ]
    
    return all(conditions)

# Iterate over all possible combinations
for num1, num2 in permutations(numbers, 2):
    for let1, let2 in permutations(letters, 2):
        if is_valid_combination(num1, num2, let1, let2):
            password = [num1, num2, let1, let2]
            print(f"<<< {password} >>>")
            break