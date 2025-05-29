from itertools import permutations

# Possible numbers and letters based on deductions
numbers = ['1', '8']
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Generate all permutations of the numbers and letters
possible_combinations = permutations(numbers + letters, 4)

# Function to check if a combination satisfies all conditions
def check_combination(combination):
    # Check against each guess
    if (combination[0] == '1' and combination[1] == '8' and combination[2] == 'H'):
        return True
    return False

# Find the correct combination
for combination in possible_combinations:
    if check_combination(combination):
        print(combination)
        break