def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(true_statements):
    # Count true and false statements
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Check each statement's validity
    results = [False] * 7
    
    # Statement 1: At least 7 statements are true
    results[0] = (num_true >= 7)
    
    # Statement 2: At most 1 statement is false
    results[1] = (num_false <= 1)
    
    # Statement 3: Exactly 3 statements are true
    results[2] = (num_true == 3)
    
    # Statement 4: Exactly 0 statements are false
    results[3] = (num_false == 0)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = ((2 in true_statements) != (3 in true_statements))
    
    # Statement 6: Number of true statements is prime
    results[5] = is_prime(num_true)
    
    # Statement 7: Number of false statements is composite
    results[6] = is_composite(num_false)
    
    # Check if the results match the given true_statements
    return set(results) == set([i in true_statements for i in range(7)])

# Try all possible combinations
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