def check_statements(combination):
    # Count true statements in the combination
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 5 statements are true
    results[0] = (true_count >= 5)
    
    # Statement 2: At most 5 statements are false
    results[1] = (false_count <= 5)
    
    # Statement 3: Exactly 5 statements are true
    results[2] = (true_count == 5)
    
    # Statement 4: Exactly 2 statements are false
    results[3] = (false_count == 2)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (results[2] != results[3]) and (results[2] or results[3])
    
    # Statement 6: The number of true statements is a prime number
    primes = {2, 3, 5, 7}
    results[5] = (true_count in primes)
    
    # Statement 7: The number of false statements is a composite number
    composites = {4, 6}  # only possible composites for false_count
    results[6] = (false_count in composites)
    
    # Check if the results match the combination
    return all(r == c for r, c in zip(results, combination))

# Generate all possible combinations and check them
solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {sol}")
        print(f"Number of true statements: {sum(sol)}")
        print(f"Number of false statements: {7-sum(sol)}")
        print()