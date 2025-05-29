def check_statements(truths):
    # Count number of true and false statements
    true_count = sum(1 for t in truths if t)
    false_count = 7 - true_count
    
    # For each statement, check if its claimed condition matches its truth value
    checks = [False] * 7
    
    # Statement 1: At least 7 of these 7 statements are true
    checks[0] = (true_count >= 7)
    
    # Statement 2: At most 2 of these 7 statements are false
    checks[1] = (false_count <= 2)
    
    # Statement 3: Exactly 2 of these 7 statements are true
    checks[2] = (true_count == 2)
    
    # Statement 4: Exactly 4 of these 7 statements are false
    checks[3] = (false_count == 4)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    checks[4] = ((truths[2] or truths[3]) and not (truths[2] and truths[3]))
    
    # Statement 6: The number of true statements is a prime number
    primes = {2, 3, 5, 7}
    checks[5] = (true_count in primes)
    
    # Statement 7: The number of false statements is a composite number
    composites = {4, 6}
    checks[6] = (false_count in composites)
    
    # A solution is valid if each statement's truth value matches its check
    return all(truths[i] == checks[i] for i in range(7))

# Test all possible combinations
solutions = []
for i in range(2**7):
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        solutions.append(truths)
        print(f"Found solution: {truths}")
        print(f"True statements: {sum(truths)}")
        print(f"False statements: {7-sum(truths)}")

print(f"Total number of solutions: {len(solutions)}")