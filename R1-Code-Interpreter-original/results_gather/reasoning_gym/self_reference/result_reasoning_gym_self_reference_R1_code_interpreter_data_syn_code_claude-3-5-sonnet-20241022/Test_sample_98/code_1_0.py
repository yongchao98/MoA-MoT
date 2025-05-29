def check_statements(truth_values):
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # For each statement, check if its truth value matches what it claims
    checks = [False] * 7
    
    # Statement 1: At least 7 are true
    checks[0] = (true_count >= 7) == truth_values[0]
    
    # Statement 2: At most 2 are false
    checks[1] = (false_count <= 2) == truth_values[1]
    
    # Statement 3: Exactly 5 are true
    checks[2] = (true_count == 5) == truth_values[2]
    
    # Statement 4: Exactly 2 are false
    checks[3] = (false_count == 2) == truth_values[3]
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    stmt3_true = truth_values[2]
    stmt4_true = truth_values[3]
    checks[4] = ((stmt3_true != stmt4_true) and (stmt3_true or stmt4_true)) == truth_values[4]
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    checks[5] = (true_count in primes) == truth_values[5]
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}
    checks[6] = (false_count in composites) == truth_values[6]
    
    # If all checks pass, this is a valid solution
    return all(checks)

# Test all possible combinations
solutions = []
for i in range(2**7):
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)
        print(f"Found solution: {truth_values}")
        print(f"True statements: {sum(truth_values)}")
        print(f"False statements: {7-sum(truth_values)}")

print(f"\nTotal number of solutions: {len(solutions)}")