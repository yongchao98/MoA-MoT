from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

# Generate all possible combinations of two numbers and two letters
all_combinations = [
    (n1, n2, l1, l2)
    for n1, n2 in permutations(possible_numbers, 2)
    for l1, l2 in permutations(possible_letters, 2)
]

# Function to check if a combination satisfies all constraints
def satisfies_constraints(combination):
    n1, n2, l1, l2 = combination
    
    # Apply constraints from each guess
    # Guess 1: 27JR
    if n1 in '27' or n2 in '27' or l1 in 'JR' or l2 in 'JR':
        return False
    
    # Guess 2: 67GP
    if n1 in '67' or n2 in '67' or (l1 == 'G' and l2 != 'P') or (l2 == 'G' and l1 != 'P'):
        return False
    
    # Guess 3: 09CO
    if (n1 == '0' and n2 == '9') or (n2 == '0' and n1 == '9') or l1 in 'CO' or l2 in 'CO':
        return False
    
    # Guess 4: 57VJ
    if n1 in '57' or n2 in '57' or l1 in 'VJ' or l2 in 'VJ':
        return False
    
    # Guess 5: 06HO
    if (n1 == '0' and n2 == '6') or (n2 == '0' and n1 == '6') or l1 in 'HO' or l2 in 'HO':
        return False
    
    # Guess 6: 67JQ
    if n1 in '67' or n2 in '67' or (l1 == 'Q' and l2 != 'J') or (l2 == 'Q' and l1 != 'J'):
        return False
    
    # Guess 7: 34TS
    if (n1 == '3' and n2 != '4') or (n2 == '3' and n1 != '4') or l1 in 'TS' or l2 in 'TS':
        return False
    
    return True

# Filter combinations that satisfy all constraints
valid_combinations = [comb for comb in all_combinations if satisfies_constraints(comb)]

# Print the valid combination(s)
for combination in valid_combinations:
    print(f"<<< {list(combination)} >>>")