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
    validities = [False] * 7
    
    # Statement 1: At least 5 are true
    validities[0] = (num_true >= 5)
    
    # Statement 2: At most 4 are false
    validities[1] = (num_false <= 4)
    
    # Statement 3: Exactly 1 is true
    validities[2] = (num_true == 1)
    
    # Statement 4: Exactly 5 are false
    validities[3] = (num_false == 5)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    validities[4] = (validities[2] != validities[3]) and (validities[2] or validities[3])
    
    # Statement 6: Number of true statements is prime
    validities[5] = is_prime(num_true)
    
    # Statement 7: Number of false statements is composite
    validities[6] = is_composite(num_false)
    
    # Check if the validities match the true_statements
    for i in range(7):
        if (i in true_statements) != validities[i]:
            return False
    return True

# Try all possible combinations
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
        print(f"True statements: {sorted(list(sol+1))}")