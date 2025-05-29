from itertools import combinations

# Total number of statements
n = 7

# Prime numbers less than or equal to 7
prime_numbers = {2, 3, 5, 7}

# Composite numbers less than or equal to 7
composite_numbers = {4, 6}

# Function to check if a combination of true statements is valid
def is_valid_combination(true_statements):
    num_true = len(true_statements)
    num_false = n - num_true
    
    # Check each statement
    if 1 in true_statements and num_true < 7:
        return False
    if 2 in true_statements and num_true < 6:
        return False
    if 3 in true_statements and num_true != 3:
        return False
    if 4 in true_statements and num_true != 7:
        return False
    if 5 in true_statements and not (3 in true_statements) ^ (4 in true_statements):
        return False
    if 6 in true_statements and num_true not in prime_numbers:
        return False
    if 7 in true_statements and num_false not in composite_numbers:
        return False
    
    return True

# Count valid combinations
valid_combinations = 0

# Check all combinations of true statements
for r in range(1, n + 1):
    for true_statements in combinations(range(1, n + 1), r):
        if is_valid_combination(true_statements):
            valid_combinations += 1

print(valid_combinations)