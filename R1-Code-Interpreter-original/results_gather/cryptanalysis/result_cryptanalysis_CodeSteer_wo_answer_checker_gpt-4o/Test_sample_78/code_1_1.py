from itertools import permutations

# Define possible numbers and letters
possible_numbers = ['3', '2']
possible_letters = ['D', 'H']

# Generate all permutations of numbers and letters
all_combinations = permutations(possible_numbers + possible_letters, 4)

# Function to validate a combination against the conditions
def is_valid_combination(combination):
    # Convert combination to a string for easier comparison
    comb_str = ''.join(combination)
    
    # Check each condition
    # Condition 1: 14HF
    if not (comb_str[0] != '1' and comb_str[1] != '4' and
            (comb_str[2] == 'H' or comb_str[3] == 'H') and
            (comb_str[2] != 'F' and comb_str[3] != 'F')):
        return False
    
    # Condition 2: 38FX
    if not ((comb_str[0] == '3' or comb_str[1] == '3') and
            (comb_str[0] != '8' and comb_str[1] != '8') and
            (comb_str[2] != 'F' and comb_str[3] != 'F') and
            (comb_str[2] != 'X' and comb_str[3] != 'X')):
        return False
    
    # Condition 3: 52YF
    if not ((comb_str[0] == '2' or comb_str[1] == '2') and
            (comb_str[0] != '5' and comb_str[1] != '5') and
            (comb_str[2] != 'Y' and comb_str[3] != 'Y') and
            (comb_str[2] != 'F' and comb_str[3] != 'F')):
        return False
    
    # Condition 4: 32QI
    if not (comb_str[0] == '2' and comb_str[1] == '3' and
            (comb_str[2] != 'Q' and comb_str[3] != 'Q') and
            (comb_str[2] != 'I' and comb_str[3] != 'I')):
        return False
    
    # Condition 5: 71DX
    if not (comb_str[0] != '7' and comb_str[1] != '1' and
            comb_str[2] == 'D' and
            (comb_str[3] != 'X')):
        return False
    
    return True

# Find the valid combination
for combination in all_combinations:
    if is_valid_combination(combination):
        print(f"<<< {list(combination)} >>>")
        break