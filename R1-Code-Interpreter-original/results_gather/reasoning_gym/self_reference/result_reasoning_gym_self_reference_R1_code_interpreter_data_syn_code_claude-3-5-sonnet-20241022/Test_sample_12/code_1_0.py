def check_statements(combination):
    # Count true and false statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 2 are true
    results[0] = (true_count >= 2)
    
    # Statement 2: At most 5 are false
    results[1] = (false_count <= 5)
    
    # Statement 3: Exactly 0 are true
    results[2] = (true_count == 0)
    
    # Statement 4: Exactly 3 are false
    results[3] = (false_count == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (combination[2] != combination[3] and (combination[2] or combination[3]))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    results[6] = (false_count in composites)
    
    # Check if the results match the combination
    return all(r == c for r, c in zip(results, combination))

# Generate and test all possible combinations
valid_solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        valid_solutions.append(combination)

print(f"Number of valid solutions: {len(valid_solutions)}")
print("\nValid solutions:")
for solution in valid_solutions:
    print(f"Solution: {solution}")
    print(f"True statements: {sum(solution)}")
    print(f"False statements: {7 - sum(solution)}")
    print("---")