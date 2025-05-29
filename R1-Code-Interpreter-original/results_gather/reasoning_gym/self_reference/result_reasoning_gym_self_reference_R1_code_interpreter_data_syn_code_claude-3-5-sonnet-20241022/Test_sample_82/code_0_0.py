def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(combination):
    # Convert to list of booleans
    statements = [bool(x) for x in combination]
    
    # Count true and false statements
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement
    checks = [
        # Statement 1: At least 3 true
        true_count >= 3,
        
        # Statement 2: At most 2 false
        false_count <= 2,
        
        # Statement 3: Exactly 6 true
        true_count == 6,
        
        # Statement 4: Exactly 1 false
        false_count == 1,
        
        # Statement 5: Statement 3 XOR Statement 4
        statements[2] != statements[3],
        
        # Statement 6: Number of true statements is prime
        is_prime(true_count),
        
        # Statement 7: Number of false statements is composite
        is_composite(false_count)
    ]
    
    # Check if each statement's truth value matches its check
    return all(check == statement for check, statement in zip(checks, statements))

# Generate all possible combinations and count valid solutions
solutions = []
for i in range(2**7):
    # Convert number to binary representation of length 7
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(sol)