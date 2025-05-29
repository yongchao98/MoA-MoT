def check_statements(true_statements):
    # Count of true statements
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Helper functions
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 4 are true
    results[0] = (num_true >= 4)
    
    # Statement 2: At most 6 are false
    results[1] = (num_false <= 6)
    
    # Statement 3: Exactly 4 are true
    results[2] = (num_true == 4)
    
    # Statement 4: Exactly 3 are false
    results[3] = (num_false == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = ((2 in true_statements) != (3 in true_statements))
    
    # Statement 6: Number of true statements is prime
    results[5] = is_prime(num_true)
    
    # Statement 7: Number of false statements is composite
    results[6] = is_composite(num_false)
    
    # Check if this combination is valid
    return set(range(7)) == true_statements.union(set([i for i in range(7) if not results[i]]))

# Generate all possible combinations
solutions = []
for i in range(128):  # 2^7 possibilities
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