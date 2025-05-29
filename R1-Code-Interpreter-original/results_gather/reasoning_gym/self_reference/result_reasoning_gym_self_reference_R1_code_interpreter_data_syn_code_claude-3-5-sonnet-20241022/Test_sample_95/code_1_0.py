def check_statements(combination):
    # Count true and false statements
    true_count = 0
    false_count = 0
    
    # First, evaluate each statement based on the combination
    # Statement 1: At least 6 are true
    s1 = true_count >= 6
    
    # Statement 2: At most 4 are false
    s2 = false_count <= 4
    
    # Statement 3: Exactly 5 are true
    s3 = true_count == 5
    
    # Statement 4: All are false
    s4 = true_count == 0
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    s5 = (s3 != s4) and (s3 or s4)
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    s6 = true_count in primes
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}
    s7 = false_count in composites
    
    # Now check if the statements match their assigned truth values
    statements = [s1, s2, s3, s4, s5, s6, s7]
    true_count = sum(1 for i in range(7) if combination[i])
    false_count = 7 - true_count
    
    # Recalculate each statement with actual counts
    s1 = true_count >= 6
    s2 = false_count <= 4
    s3 = true_count == 5
    s4 = true_count == 0
    s5 = (s3 != s4) and (s3 or s4)
    s6 = true_count in primes
    s7 = false_count in composites
    
    statements = [s1, s2, s3, s4, s5, s6, s7]
    
    # Check if the combination matches what each statement claims
    return all(statements[i] == combination[i] for i in range(7))

# Test all possible combinations
solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("Solutions found:")
    for sol in solutions:
        print(sol)