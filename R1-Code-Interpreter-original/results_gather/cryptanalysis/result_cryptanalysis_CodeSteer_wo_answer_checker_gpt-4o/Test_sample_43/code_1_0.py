from itertools import permutations

# Possible numbers and letters based on the conditions
possible_numbers = ['3', '5', '8']
possible_letters = ['P', 'T']

# Generate all permutations of the possible numbers and letters
number_permutations = permutations(possible_numbers, 2)
letter_permutations = permutations(possible_letters, 2)

# Function to validate a combination against the conditions
def is_valid_combination(numbers, letters):
    # Check each condition
    # Condition 1: 61OJ - both numbers and letters are incorrect
    if numbers[0] in ['6', '1'] or numbers[1] in ['6', '1']:
        return False
    if letters[0] in ['O', 'J'] or letters[1] in ['O', 'J']:
        return False
    
    # Condition 2: 98QN - one number correct and in the correct position
    if numbers[0] == '9' or numbers[1] == '9':
        return False
    if letters[0] in ['Q', 'N'] or letters[1] in ['Q', 'N']:
        return False
    
    # Condition 3: 58FC - one number correct and in the correct position
    if numbers[0] != '5' or numbers[1] == '5':
        return False
    if letters[0] in ['F', 'C'] or letters[1] in ['F', 'C']:
        return False
    
    # Condition 4: 72YK - both numbers and letters are incorrect
    if numbers[0] in ['7', '2'] or numbers[1] in ['7', '2']:
        return False
    if letters[0] in ['Y', 'K'] or letters[1] in ['Y', 'K']:
        return False
    
    # Condition 5: 02ER - both numbers and letters are incorrect
    if numbers[0] in ['0', '2'] or numbers[1] in ['0', '2']:
        return False
    if letters[0] in ['E', 'R'] or letters[1] in ['E', 'R']:
        return False
    
    # Condition 6: 35QW - one number correct and in the correct position
    if numbers[0] != '5' or numbers[1] == '5':
        return False
    if letters[0] in ['Q', 'W'] or letters[1] in ['Q', 'W']:
        return False
    
    # Condition 7: 71KJ - both numbers and letters are incorrect
    if numbers[0] in ['7', '1'] or numbers[1] in ['7', '1']:
        return False
    if letters[0] in ['K', 'J'] or letters[1] in ['K', 'J']:
        return False
    
    # Condition 8: 79PX - one letter correct and in the correct position
    if numbers[0] in ['7', '9'] or numbers[1] in ['7', '9']:
        return False
    if letters[0] != 'P' or letters[1] == 'P':
        return False
    
    # Condition 9: 86TV - one number correct but in the wrong position
    if numbers[0] == '8' or numbers[1] != '8':
        return False
    if letters[0] != 'T' or letters[1] == 'T':
        return False
    
    return True

# Iterate through all permutations and find the valid combination
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        if is_valid_combination(num_perm, let_perm):
            password = list(num_perm) + list(let_perm)
            print(f"<<< {password} >>>")
            break