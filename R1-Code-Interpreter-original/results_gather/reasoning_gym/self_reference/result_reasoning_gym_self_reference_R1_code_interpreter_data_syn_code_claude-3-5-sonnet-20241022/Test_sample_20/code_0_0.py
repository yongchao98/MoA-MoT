def check_statements(true_statements):
    # Count of true and false statements
    true_count = len(true_statements)
    false_count = 7 - true_count
    
    # Create a list of booleans indicating which statements are true
    is_true = [i in true_statements for i in range(7)]
    
    # Check each statement
    valid = [False] * 7
    
    # Statement 1: At least 3 are true
    valid[0] = (true_count >= 3)
    
    # Statement 2: At most 3 are false
    valid[1] = (false_count <= 3)
    
    # Statement 3: Exactly 7 are true
    valid[2] = (true_count == 7)
    
    # Statement 4: Exactly 2 are false
    valid[3] = (false_count == 2)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid[4] = (is_true[2] != is_true[3]) and (is_true[2] or is_true[3])
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    valid[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    valid[6] = (false_count in composites)
    
    # Check if the combination is consistent
    return all(x == is_true[i] for i, x in enumerate(valid))

# Generate all possible combinations and check them
solutions = []
for i in range(128):  # 2^7 = 128 combinations
    true_statements = set()
    for j in range(7):
        if (i >> j) & 1:
            true_statements.add(j)
    if check_statements(true_statements):
        solutions.append(true_statements)

print(len(solutions))