from itertools import permutations

# Possible numbers and letters
numbers = ['0', '1']
letters = ['A', 'B', 'C', 'D', 'H', 'J']

# Generate all permutations of the numbers and letters
possible_combinations = permutations(numbers + letters, 4)

# Function to check if a combination satisfies all conditions
def check_combination(combination):
    num1, num2, let1, let2 = combination
    
    # Check each condition
    # Condition 1: 15PU
    if not ((num1 == '1' or num2 == '1') and (num1 != '5' and num2 != '5') and let1 < 'P' and let2 < 'U'):
        return False
    
    # Condition 2: 43EH
    if not (num1 < '4' and num2 < '3' and (let1 == 'E' or let2 == 'E') and (let1 != 'H' and let2 != 'H')):
        return False
    
    # Condition 3: 94HG
    if not (num1 < '9' and num2 < '4' and (let1 == 'H' or let2 == 'H') and (let1 != 'G' and let2 != 'G')):
        return False
    
    # Condition 4: 36NW
    if not (num1 < '3' and num2 < '6' and let1 < 'N' and let2 < 'W'):
        return False
    
    # Condition 5: 39FU
    if not (num1 < '3' and num2 < '9' and let1 != 'F' and let2 != 'U'):
        return False
    
    # Condition 6: 27JZ
    if not (num1 < '2' and num2 < '7' and (let1 == 'J' or let2 == 'J') and (let1 != 'Z' and let2 != 'Z')):
        return False
    
    return True

# Find the correct combination
for combination in possible_combinations:
    if check_combination(combination):
        print(combination)
        break