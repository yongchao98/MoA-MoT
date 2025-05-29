from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Generate all possible combinations of two numbers and two letters
possible_combinations = [''.join(p + q) for p in permutations(numbers, 2) for q in permutations(letters, 2)]

# Function to check if a combination satisfies the constraints
def satisfies_constraints(combination):
    num1, num2, let1, let2 = combination

    # Guess 74JY
    if not ((num1 == '7' or num2 == '7') and (num1 != '4' and num2 != '4')): return False
    if not ((let1 == 'J' or let2 == 'J') and (let1 != 'Y' and let2 != 'Y')): return False
    if let1 <= 'J' or let2 <= 'J': return False

    # Guess 93ZN
    if num1 in '93' or num2 in '93' or let1 in 'ZN' or let2 in 'ZN': return False

    # Guess 26MU
    if num1 in '0126' or num2 in '0126' or let1 in 'MU' or let2 in 'MU': return False

    # Guess 57FS
    if not ((num1 == '5' or num2 == '5') and (num1 != '7' and num2 != '7')): return False
    if not ((let1 == 'S' or let2 == 'S') and (let1 != 'F' and let2 != 'F')): return False
    if let1 <= 'S' or let2 <= 'S': return False

    return True

# Filter combinations that satisfy all constraints
valid_combinations = [comb for comb in possible_combinations if satisfies_constraints(comb)]

# Output the valid combination
if valid_combinations:
    password = list(valid_combinations[0])
    print(f"<<< {password} >>>")
else:
    print("No valid combination found.")