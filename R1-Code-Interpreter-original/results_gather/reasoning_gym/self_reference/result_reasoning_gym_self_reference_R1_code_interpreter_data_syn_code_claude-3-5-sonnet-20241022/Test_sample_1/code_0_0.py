def check_statements(combination):
    # Count true and false statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement's validity
    results = [False] * 7
    
    # Statement 1: At least 6 statements are true
    results[0] = (true_count >= 6) == combination[0]
    
    # Statement 2: At most 7 statements are false
    results[1] = (false_count <= 7) == combination[1]
    
    # Statement 3: Exactly 5 statements are true
    results[2] = (true_count == 5) == combination[2]
    
    # Statement 4: Exactly 2 statements are false
    results[3] = (false_count == 2) == combination[3]
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    stmt3_true = combination[2]
    stmt4_true = combination[3]
    results[4] = ((stmt3_true != stmt4_true) and (stmt3_true or stmt4_true)) == combination[4]
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (true_count in primes) == combination[5]
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    results[6] = (false_count in composites) == combination[6]
    
    # Check if all statements are consistent
    return all(results)

# Generate and check all possible combinations
valid_solutions = 0
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        valid_solutions += 1

print(valid_solutions)