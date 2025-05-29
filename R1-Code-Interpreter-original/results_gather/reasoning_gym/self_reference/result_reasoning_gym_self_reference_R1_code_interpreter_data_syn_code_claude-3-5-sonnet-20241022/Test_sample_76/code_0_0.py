def check_statements(combination):
    # Count true statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 2 are true
    results[0] = (true_count >= 2)
    
    # Statement 2: At most 2 are false
    results[1] = (false_count <= 2)
    
    # Statement 3: Exactly 4 are true
    results[2] = (true_count == 4)
    
    # Statement 4: Exactly 1 is false
    results[3] = (false_count == 1)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (combination[2] != combination[3]) and (combination[2] or combination[3])
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composite = {4, 6}  # only possible composite numbers for false statements
    results[6] = (false_count in composite)
    
    # Check if the combination matches the results
    return all(combination[i] == results[i] for i in range(7))

# Generate and test all possible combinations
solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    combination = [(i >> j) & 1 == 1 for j in range(7)]
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