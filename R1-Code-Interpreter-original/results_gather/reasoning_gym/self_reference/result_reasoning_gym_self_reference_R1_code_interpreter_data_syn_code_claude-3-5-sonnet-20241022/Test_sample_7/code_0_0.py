def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(statements):
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement's validity
    results = [
        true_count >= 6,                    # Statement 1
        false_count <= 0,                   # Statement 2
        true_count == 2,                    # Statement 3
        false_count == 1,                   # Statement 4
        (statements[2] != statements[3]),   # Statement 5
        is_prime(true_count),              # Statement 6
        is_composite(false_count)          # Statement 7
    ]
    
    # Check if the results match the input statements
    return results == list(statements)

# Try all possible combinations
solutions = []
for i in range(2**7):
    statements = [bool((i >> j) & 1) for j in range(7)]
    if check_statements(statements):
        solutions.append(statements)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Statements: {[int(x) for x in sol]}")
        print(f"True count: {sum(sol)}")
        print(f"False count: {7-sum(sol)}")
        print()