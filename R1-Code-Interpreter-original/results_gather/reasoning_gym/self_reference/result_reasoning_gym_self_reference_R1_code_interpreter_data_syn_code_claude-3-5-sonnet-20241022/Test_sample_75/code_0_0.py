def check_statements(truths):
    # Count true and false statements
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # Check each statement's validity
    valid = [False] * 7
    
    # Statement 1: At least 7 statements are true
    valid[0] = (true_count >= 7) == truths[0]
    
    # Statement 2: At most 3 statements are false
    valid[1] = (false_count <= 3) == truths[1]
    
    # Statement 3: Exactly 4 statements are true
    valid[2] = (true_count == 4) == truths[2]
    
    # Statement 4: Exactly 7 statements are false
    valid[3] = (false_count == 7) == truths[3]
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid[4] = ((truths[2] != truths[3]) and (truths[2] or truths[3])) == truths[4]
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    valid[5] = (true_count in primes) == truths[5]
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    valid[6] = (false_count in composites) == truths[6]
    
    return all(valid)

# Try all possible combinations
solutions = []
for i in range(2**7):
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        solutions.append(truths)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {sol}")
        print(f"Number of true statements: {sum(sol)}")
        print(f"Number of false statements: {7-sum(sol)}")
        print()