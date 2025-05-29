def check_statements(truths):
    # Count true and false statements
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 6 statements are true
    results[0] = (true_count >= 6)
    
    # Statement 2: At most 4 statements are false
    results[1] = (false_count <= 4)
    
    # Statement 3: Exactly 5 statements are true
    results[2] = (true_count == 5)
    
    # Statement 4: Exactly 2 statements are false
    results[3] = (false_count == 2)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (truths[2] != truths[3]) and (truths[2] or truths[3])
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for false_count
    results[6] = (false_count in composites)
    
    # Check if this combination is valid
    return all(truths[i] == results[i] for i in range(7))

# Generate all possible combinations and check them
solutions = []
for i in range(2**7):
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        solutions.append(truths)

print(len(solutions))