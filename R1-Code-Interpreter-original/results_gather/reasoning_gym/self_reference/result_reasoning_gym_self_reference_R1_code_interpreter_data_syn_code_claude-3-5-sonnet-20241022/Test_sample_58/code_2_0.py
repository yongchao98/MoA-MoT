def check_statements(combination):
    # Count true and false statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # First verify if the statements match their claimed truth values
    
    # Statement 1: At least 2 are true
    if combination[0] and true_count < 2: return False
    if not combination[0] and true_count >= 2: return False
    
    # Statement 2: At most 4 are false
    if combination[1] and false_count > 4: return False
    if not combination[1] and false_count <= 4: return False
    
    # Statement 3: Exactly 5 are true
    if combination[2] and true_count != 5: return False
    if not combination[2] and true_count == 5: return False
    
    # Statement 4: Exactly 2 are false
    if combination[3] and false_count != 2: return False
    if not combination[3] and false_count == 2: return False
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    stmt5_should_be = (combination[2] != combination[3] and (combination[2] or combination[3]))
    if combination[4] != stmt5_should_be: return False
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    if combination[5] and true_count not in primes: return False
    if not combination[5] and true_count in primes: return False
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}
    if combination[6] and false_count not in composites: return False
    if not combination[6] and false_count in composites: return False
    
    return True

solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)
        print(f"Solution found: {combination}")
        print(f"Number of true statements: {sum(combination)}")
        print(f"Number of false statements: {7-sum(combination)}")
        print("---")

print(f"Total number of solutions: {len(solutions)}")