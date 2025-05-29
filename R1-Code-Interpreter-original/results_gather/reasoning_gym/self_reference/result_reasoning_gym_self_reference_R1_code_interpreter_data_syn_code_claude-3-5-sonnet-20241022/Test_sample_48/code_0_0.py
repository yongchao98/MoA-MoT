def check_statements(truths):
    # Count number of true and false statements
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # Check each statement's consistency
    consistent = [False] * 7
    
    # Statement 1: At least 7 statements are true
    consistent[0] = (truths[0] == (true_count >= 7))
    
    # Statement 2: At most 2 statements are false
    consistent[1] = (truths[1] == (false_count <= 2))
    
    # Statement 3: Exactly 2 statements are true
    consistent[2] = (truths[2] == (true_count == 2))
    
    # Statement 4: Exactly 4 statements are false
    consistent[3] = (truths[3] == (false_count == 4))
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    consistent[4] = (truths[4] == (truths[2] != truths[3] and (truths[2] or truths[3])))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    consistent[5] = (truths[5] == (true_count in primes))
    
    # Statement 6: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    consistent[6] = (truths[6] == (false_count in composites))
    
    # Check if all statements are consistent
    return all(consistent[i] == truths[i] for i in range(7))

# Generate all possible combinations and check them
solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        solutions.append(truths)

print(len(solutions))