def check_statements(truths):
    # Count true and false statements
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # For each statement, check if its truth value matches what it claims
    checks = [False] * 7
    
    # Statement 1: At least 7 of these 7 statements are true
    checks[0] = (true_count >= 7)
    
    # Statement 2: At most 5 of these 7 statements are false
    checks[1] = (false_count <= 5)
    
    # Statement 3: Exactly 4 of these 7 statements are true
    checks[2] = (true_count == 4)
    
    # Statement 4: Exactly 3 of these 7 statements are false
    checks[3] = (false_count == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    checks[4] = (truths[2] != truths[3] and (truths[2] or truths[3]))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    checks[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}
    checks[6] = (false_count in composites)
    
    # A combination is valid if each statement's truth value matches what it claims
    return all(checks[i] == truths[i] for i in range(7))

# Try all possible combinations
solutions = []
for i in range(2**7):
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        solutions.append(truths)
        print(f"Found solution: {truths}")
        print(f"Number of true statements: {sum(truths)}")
        print(f"Number of false statements: {7-sum(truths)}")

print(f"\nTotal number of solutions: {len(solutions)}")