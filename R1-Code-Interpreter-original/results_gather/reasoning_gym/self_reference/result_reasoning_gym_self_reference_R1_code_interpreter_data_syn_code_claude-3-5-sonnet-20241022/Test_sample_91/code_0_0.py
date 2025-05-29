def check_statements(statements):
    # Count true statements
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement's validity
    valid = [False] * 7
    
    # Statement 1: At least 4 are true
    valid[0] = (true_count >= 4) == statements[0]
    
    # Statement 2: At most 0 are false
    valid[1] = (false_count <= 0) == statements[1]
    
    # Statement 3: Exactly 7 are true
    valid[2] = (true_count == 7) == statements[2]
    
    # Statement 4: Exactly 3 are false
    valid[3] = (false_count == 3) == statements[3]
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid[4] = ((statements[2] != statements[3]) and (statements[2] or statements[3])) == statements[4]
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    valid[5] = (true_count in primes) == statements[5]
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    valid[6] = (false_count in composites) == statements[6]
    
    # Check if all statements are consistent
    return all(valid)

# Try all possible combinations
solutions = []
for i in range(2**7):
    statements = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(statements):
        solutions.append(statements)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Statements: {[1 if x else 0 for x in sol]}")
        print(f"True count: {sum(sol)}")
        print("---")