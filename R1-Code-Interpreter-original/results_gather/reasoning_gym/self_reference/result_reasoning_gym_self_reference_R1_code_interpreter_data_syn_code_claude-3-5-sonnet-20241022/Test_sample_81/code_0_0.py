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
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Create array of statement truth values
    statements = [i in true_statements for i in range(7)]
    
    # Check each statement's claim
    checks = [
        num_true >= 7,                    # Statement 1
        num_false <= 5,                   # Statement 2
        num_true == 3,                    # Statement 3
        num_false == 0,                   # Statement 4
        (statements[2] != statements[3]), # Statement 5
        is_prime(num_true),              # Statement 6
        is_composite(num_false)          # Statement 7
    ]
    
    # Statement is consistent if its truth value matches its claim
    for i in range(7):
        if statements[i] != checks[i]:
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