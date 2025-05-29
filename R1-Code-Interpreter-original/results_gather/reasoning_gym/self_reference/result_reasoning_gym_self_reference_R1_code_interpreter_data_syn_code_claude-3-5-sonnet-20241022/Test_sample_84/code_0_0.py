def check_statements(true_statements):
    # Count of true statements
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 3 are true
    results[0] = (num_true >= 3)
    
    # Statement 2: At most 3 are false
    results[1] = (num_false <= 3)
    
    # Statement 3: Exactly 4 are true
    results[2] = (num_true == 4)
    
    # Statement 4: Exactly 3 are false
    results[3] = (num_false == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    stmt3_true = (2 in true_statements)
    stmt4_true = (3 in true_statements)
    results[4] = (stmt3_true != stmt4_true)
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (num_true in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for false statements (1-7)
    results[6] = (num_false in composites)
    
    # Check if the results match the true_statements
    return all(i in true_statements for i, r in enumerate(results) if r) and \
           all(i not in true_statements for i, r in enumerate(results) if not r)

# Generate all possible combinations
solutions = []
for i in range(128):  # 2^7 = 128 combinations
    true_statements = set()
    for j in range(7):
        if (i >> j) & 1:
            true_statements.add(j)
    if check_statements(true_statements):
        solutions.append(true_statements)

print(f"Number of solutions: {len(solutions)}")
if len(solutions) > 0:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"True statements: {sorted(list(sol))}")