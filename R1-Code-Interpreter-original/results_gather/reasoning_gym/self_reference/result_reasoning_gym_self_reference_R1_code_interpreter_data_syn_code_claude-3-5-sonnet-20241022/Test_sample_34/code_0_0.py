def check_statements(statements):
    # Count number of true statements
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement's consistency
    consistent = [False] * 7
    
    # Statement 1: At least 4 are true
    consistent[0] = (statements[0] == (true_count >= 4))
    
    # Statement 2: At most 3 are false
    consistent[1] = (statements[1] == (false_count <= 3))
    
    # Statement 3: Exactly 2 are true
    consistent[2] = (statements[2] == (true_count == 2))
    
    # Statement 4: All are false
    consistent[3] = (statements[3] == (true_count == 0))
    
    # Statement 5: Statement 3 XOR Statement 4 is true
    consistent[4] = (statements[4] == (statements[2] != statements[3] and (statements[2] or statements[3])))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    consistent[5] = (statements[5] == (true_count in primes))
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # composites from 1-7
    consistent[6] = (statements[6] == (false_count in composites))
    
    # Check if all statements are consistent
    return all(consistent)

# Generate all possible combinations and check each
solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    statements = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(statements):
        solutions.append(statements)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"True statements: {[i+1 for i in range(7) if sol[i]]}")