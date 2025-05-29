def check_statements(truth_values):
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # Check each statement's consistency
    consistent = [False] * 7
    
    # Statement 1: At least 7 statements are true
    consistent[0] = (truth_values[0] == (true_count >= 7))
    
    # Statement 2: At most 2 statements are false
    consistent[1] = (truth_values[1] == (false_count <= 2))
    
    # Statement 3: Exactly 5 statements are true
    consistent[2] = (truth_values[2] == (true_count == 5))
    
    # Statement 4: Exactly 2 statements are false
    consistent[3] = (truth_values[3] == (false_count == 2))
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    consistent[4] = (truth_values[4] == (truth_values[2] != truth_values[3] and (truth_values[2] or truth_values[3])))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    consistent[5] = (truth_values[5] == (true_count in primes))
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    consistent[6] = (truth_values[6] == (false_count in composites))
    
    # Check if all statements are consistent
    return all(consistent[i] == truth_values[i] for i in range(7))

# Generate all possible combinations and count valid solutions
solutions = []
for i in range(2**7):
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)

print(len(solutions))