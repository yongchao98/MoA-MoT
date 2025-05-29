def check_statements(truth_values):
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # Check each statement's validity
    valid = [False] * 7
    
    # Statement 1: At least 6 statements are true
    valid[0] = (true_count >= 6) == truth_values[0]
    
    # Statement 2: At most 3 statements are false
    valid[1] = (false_count <= 3) == truth_values[1]
    
    # Statement 3: Exactly 4 statements are true
    valid[2] = (true_count == 4) == truth_values[2]
    
    # Statement 4: Exactly 3 statements are false
    valid[3] = (false_count == 3) == truth_values[3]
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid[4] = ((truth_values[2] != truth_values[3]) and (truth_values[2] or truth_values[3])) == truth_values[4]
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    valid[5] = (true_count in primes) == truth_values[5]
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    valid[6] = (false_count in composites) == truth_values[6]
    
    # Check if all statements are consistent
    return all(valid)

# Try all possible combinations
solutions = []
for i in range(2**7):
    # Convert number to binary representation of truth values
    truth_values = [bool((i >> j) & 1) for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)

print(f"Number of solutions: {len(solutions)}")