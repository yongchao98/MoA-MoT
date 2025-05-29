def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def validate_solution(statements):
    # Convert to list of booleans
    is_true = [i in statements for i in range(7)]
    num_true = sum(is_true)
    num_false = 7 - num_true
    
    # Check each statement individually
    checks = [
        # Statement 1: At least 7 are true
        (num_true >= 7) == is_true[0],
        
        # Statement 2: At most 1 is false
        (num_false <= 1) == is_true[1],
        
        # Statement 3: Exactly 3 are true
        (num_true == 3) == is_true[2],
        
        # Statement 4: Exactly 0 are false
        (num_false == 0) == is_true[3],
        
        # Statement 5: Either Statement 3 or Statement 4 is true, but not both
        (is_true[2] != is_true[3]) == is_true[4],
        
        # Statement 6: Number of true statements is prime
        is_prime(num_true) == is_true[5],
        
        # Statement 7: Number of false statements is composite
        is_composite(num_false) == is_true[6]
    ]
    
    return all(checks)

# Try all possible combinations
solutions = []
for i in range(2**7):
    true_statements = set()
    for j in range(7):
        if (i >> j) & 1:
            true_statements.add(j)
    if validate_solution(true_statements):
        solutions.append(true_statements)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"True statements: {sorted(list(sol))}, "
              f"Number of true statements: {len(sol)}, "
              f"Number of false statements: {7-len(sol)}")