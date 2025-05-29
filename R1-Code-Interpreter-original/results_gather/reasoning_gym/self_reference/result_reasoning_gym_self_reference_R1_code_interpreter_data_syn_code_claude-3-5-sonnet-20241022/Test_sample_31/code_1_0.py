def check_statements(truth_values):
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # Initialize what each statement claims
    claims = [False] * 7
    
    # Statement 1: At least 7 of these 7 statements are true
    claims[0] = true_count >= 7
    
    # Statement 2: At most 3 of these 7 statements are false
    claims[1] = false_count <= 3
    
    # Statement 3: Exactly 4 of these 7 statements are true
    claims[2] = true_count == 4
    
    # Statement 4: Exactly 3 of these 7 statements are false
    claims[3] = false_count == 3
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    claims[4] = (truth_values[2] != truth_values[3]) and (truth_values[2] or truth_values[3])
    
    # Statement 6: The number of true statements is a prime number
    primes = {2, 3, 5, 7}
    claims[5] = true_count in primes
    
    # Statement 7: The number of false statements is a composite number
    composites = {4, 6}
    claims[6] = false_count in composites
    
    # A statement is consistent if its truth value matches what it claims
    consistent = all(truth_values[i] == claims[i] for i in range(7))
    
    if consistent:
        print(f"Found solution: {truth_values}")
        print(f"True count: {true_count}")
        print(f"Claims: {claims}")
    
    return consistent

# Test all possible combinations
solutions = []
for i in range(2**7):
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)

print(f"Number of solutions: {len(solutions)}")