from itertools import permutations

# Possible numbers and letters based on the feedback
possible_numbers = ['8', '9']
possible_letters = ['J', 'A']

# Generate all permutations of the numbers and letters
all_combinations = permutations(possible_numbers + possible_letters, 4)

# Function to check if a combination satisfies all constraints
def is_valid_combination(combination):
    # Convert combination to a string for easier comparison
    comb_str = ''.join(combination)
    
    # Check each guess against the feedback
    # Guess 1: 54VU
    if not (comb_str[0] > '5' and comb_str[1] > '4' and ('V' not in comb_str) and ('U' in comb_str)):
        return False
    
    # Guess 2: 19LA
    if not (comb_str[0] != '1' and comb_str[1] != '9' and ('L' not in comb_str) and ('A' in comb_str)):
        return False
    
    # Guess 3: 38JV
    if not (('3' not in comb_str) and ('8' in comb_str) and ('J' in comb_str) and ('V' not in comb_str)):
        return False
    
    # Guess 4: 14EC
    if not (comb_str[0] > '1' and comb_str[1] > '4' and ('E' not in comb_str) and ('C' not in comb_str)):
        return False
    
    # Guess 5: 47QU
    if not (comb_str[0] == '9' and comb_str[1] > '7' and ('Q' not in comb_str) and ('U' not in comb_str)):
        return False
    
    return True

# Find the valid combination
for combination in all_combinations:
    if is_valid_combination(combination):
        print(f"<<< {list(combination)} >>>")
        break