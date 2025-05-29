def check_statements(combination):
    # Count true statements in the combination
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 4 are true
    results[0] = (true_count >= 4)
    
    # Statement 2: At most 4 are false
    results[1] = (false_count <= 4)
    
    # Statement 3: Exactly 4 are true
    results[2] = (true_count == 4)
    
    # Statement 4: Exactly 3 are false
    results[3] = (false_count == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (combination[2] != combination[3] and (combination[2] or combination[3]))
    
    # Statement 6: Number of true statements is prime
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    results[5] = is_prime(true_count)
    
    # Statement 7: Number of false statements is composite
    def is_composite(n):
        return n > 1 and not is_prime(n)
    results[6] = is_composite(false_count)
    
    # Check if the results match the combination
    return all(r == c for r, c in zip(results, combination))

# Generate all possible combinations and check them
valid_solutions = []
for i in range(2**7):
    # Convert number to binary combination
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        valid_solutions.append(combination)

print(f"Number of valid solutions: {len(valid_solutions)}")
print("\nValid solutions:")
for solution in valid_solutions:
    print(f"Solution: {solution}")
    print(f"True statements: {sum(solution)}")
    print(f"False statements: {7-sum(solution)}")
    print("---")