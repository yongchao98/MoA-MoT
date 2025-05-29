from itertools import combinations

# Total number of statements
n = 7

# Possible prime numbers for true statements
prime_true_counts = {2, 3, 5, 7}

# Possible composite numbers for false statements
composite_false_counts = {4, 6}

# Function to check if a combination of true statements is valid
def is_valid_combination(true_statements):
    # Count the number of true and false statements
    true_count = len(true_statements)
    false_count = n - true_count
    
    # Check conditions
    if true_count in prime_true_counts and false_count in composite_false_counts:
        # Check each statement
        s1 = true_count >= 2
        s2 = false_count <= 0
        s3 = true_count == 1
        s4 = false_count == 7
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = true_count in prime_true_counts
        s7 = false_count in composite_false_counts
        
        # Return if all conditions are satisfied
        return s1 and not s2 and not s3 and not s4 and s5 and s6 and s7
    return False

# Generate all combinations of true statements
possible_solutions = 0
for r in range(1, n + 1):
    for true_statements in combinations(range(1, n + 1), r):
        if is_valid_combination(true_statements):
            possible_solutions += 1

print(possible_solutions)