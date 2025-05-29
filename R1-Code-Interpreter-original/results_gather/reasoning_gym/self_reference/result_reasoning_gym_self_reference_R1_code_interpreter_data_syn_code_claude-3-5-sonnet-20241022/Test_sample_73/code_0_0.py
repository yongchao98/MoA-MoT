def check_statements(true_statements):
    # Convert set of true statement indices to boolean list
    statements = [False] * 7
    for i in true_statements:
        statements[i] = True
    
    num_true = sum(statements)
    num_false = 7 - num_true
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 5 are true
    results[0] = (num_true >= 5)
    
    # Statement 2: At most 6 are false
    results[1] = (num_false <= 6)
    
    # Statement 3: Exactly 5 are true
    results[2] = (num_true == 5)
    
    # Statement 4: Exactly 2 are false
    results[3] = (num_false == 2)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (statements[2] != statements[3]) and (statements[2] or statements[3])
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (num_true in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    results[6] = (num_false in composites)
    
    # Check if the results match the assumed true/false pattern
    return all(i in true_statements for i in range(7) if results[i]) and \
           all(i not in true_statements for i in range(7) if not results[i])

# Generate all possible combinations
solutions = []
for i in range(2**7):
    true_statements = set()
    for j in range(7):
        if (i >> j) & 1:
            true_statements.add(j)
    if check_statements(true_statements):
        solutions.append(true_statements)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"True statements: {sorted(list(sol))}")