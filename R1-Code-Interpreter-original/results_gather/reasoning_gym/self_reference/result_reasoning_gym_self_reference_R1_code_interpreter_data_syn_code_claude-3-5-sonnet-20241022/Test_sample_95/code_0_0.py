def check_statements(combination):
    # Count true and false statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement's validity
    valid = [False] * 7
    
    # Statement 1: At least 6 statements are true
    valid[0] = (combination[0] == (true_count >= 6))
    
    # Statement 2: At most 4 statements are false
    valid[1] = (combination[1] == (false_count <= 4))
    
    # Statement 3: Exactly 5 statements are true
    valid[2] = (combination[2] == (true_count == 5))
    
    # Statement 4: All statements are false
    valid[3] = (combination[3] == (true_count == 0))
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid[4] = (combination[4] == (combination[2] ^ combination[3]))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    valid[5] = (combination[5] == (true_count in primes))
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    valid[6] = (combination[6] == (false_count in composites))
    
    # Check if all statements are consistent
    return all(valid[i] == combination[i] for i in range(7))

# Generate all possible combinations and count valid solutions
solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")