from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
def check_combination(combo):
    num1, num2, let1, let2 = combo

    # Guess 56QS
    if not ((num1 == '5' or num2 == '5') and (num1 != '6' and num2 != '6') and let1 < 'Q' and let2 < 'S'):
        return False

    # Guess 47KB
    if num1 in '47' or num2 in '47' or let1 in 'KB' or let2 in 'KB':
        return False

    # Guess 83CN
    if not ((num1 == '8' or num2 == '8') and (num1 != '3' and num2 != '3') and (let1 == 'C' or let2 == 'C') and let1 < 'N' and let2 < 'N'):
        return False

    # Guess 35JX
    if not ((num1 == '5' or num2 == '5') and (num1 != '3' and num2 != '3') and let1 < 'J' and let2 < 'X'):
        return False

    # Guess 95FG
    if not ((num1 == '5' or num2 == '5') and (num1 != '9' and num2 != '9') and (let1 == 'F' or let2 == 'F') and let1 != 'G' and let2 != 'G'):
        return False

    return True

# Generate all permutations of two numbers and two letters
for num_combo in permutations(numbers, 2):
    for let_combo in permutations(letters, 2):
        combo = num_combo + let_combo
        if check_combination(combo):
            print(f"<<< {list(combo)} >>>")
            break