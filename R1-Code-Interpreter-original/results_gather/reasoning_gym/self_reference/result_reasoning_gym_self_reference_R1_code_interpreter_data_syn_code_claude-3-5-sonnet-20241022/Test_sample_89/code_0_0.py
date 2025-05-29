def check_statements(truth_values):
    # Count true statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # Check each statement's consistency
    results = [False] * 7
    
    # Statement 1: At least 3 are true
    results[0] = (true_count >= 3)
    
    # Statement 2: At most 4 are false
    results[1] = (false_count <= 4)
    
    # Statement 3: Exactly 7 are true
    results[2] = (true_count == 7)
    
    # Statement 4: Exactly 3 are false
    results[3] = (false_count == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (truth_values[2] != truth_values[3]) and (truth_values[2] or truth_values[3])
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for false_count
    results[6] = (false_count in composites)
    
    # Check if the results match the assumed truth values
    return all(r == t for r, t in zip(results, truth_values))

# Generate and test all possible combinations
solutions = []
for i in range(2**7):
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {sol}")
        print(f"Number of true statements: {sum(sol)}")
        print("---")