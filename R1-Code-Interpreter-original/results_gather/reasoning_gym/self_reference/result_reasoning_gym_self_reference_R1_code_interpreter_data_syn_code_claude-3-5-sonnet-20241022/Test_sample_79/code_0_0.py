def check_statements(statements):
    # Count true and false statements
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement's validity
    valid = [False] * 7
    
    # Statement 1: At least 6 statements are true
    valid[0] = (statements[0] == (true_count >= 6))
    
    # Statement 2: At most 1 statement is false
    valid[1] = (statements[1] == (false_count <= 1))
    
    # Statement 3: Exactly 2 statements are true
    valid[2] = (statements[2] == (true_count == 2))
    
    # Statement 4: Exactly 3 statements are false
    valid[3] = (statements[3] == (false_count == 3))
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid[4] = (statements[4] == (statements[2] ^ statements[3]))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    valid[5] = (statements[5] == (true_count in primes))
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    valid[6] = (statements[6] == (false_count in composites))
    
    # Check if all statements are consistent
    return all(valid)

# Try all possible combinations
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
        print(f"Number of true statements: {sum(sol)}")
        print("---")