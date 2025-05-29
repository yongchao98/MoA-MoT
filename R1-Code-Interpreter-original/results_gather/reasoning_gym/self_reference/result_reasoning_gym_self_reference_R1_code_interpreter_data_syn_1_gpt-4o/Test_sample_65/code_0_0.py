from itertools import combinations

# Total number of statements
n_statements = 7

# Possible prime numbers for true statements
prime_true_counts = {2, 3, 5, 7}

# Possible composite numbers for false statements
composite_false_counts = {4, 6}

# Function to check if a combination of true statements is valid
def is_valid_combination(true_statements):
    n_true = len(true_statements)
    n_false = n_statements - n_true
    
    # Check each statement
    s1 = n_true >= 5
    s2 = n_false == 0
    s3 = n_true == 5
    s4 = n_false == 5
    s5 = (s3 or s4) and not (s3 and s4)
    s6 = n_true in prime_true_counts
    s7 = n_false in composite_false_counts
    
    # Create a list of truth values for each statement
    truth_values = [s1, s2, s3, s4, s5, s6, s7]
    
    # Check if the combination satisfies the truth values
    return all(truth_values[i] for i in true_statements)

# Count the number of valid solutions
valid_solutions = 0

# Check all combinations of true statements
for r in range(1, n_statements + 1):
    for true_statements in combinations(range(n_statements), r):
        if is_valid_combination(true_statements):
            valid_solutions += 1

print(valid_solutions)