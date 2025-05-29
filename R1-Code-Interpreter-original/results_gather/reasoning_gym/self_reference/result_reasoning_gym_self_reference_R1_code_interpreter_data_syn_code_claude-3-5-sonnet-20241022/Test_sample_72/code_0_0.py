from itertools import product

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
    # Convert combination to list of booleans
    statements = list(combination)
    
    # Count true and false statements
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement
    checks = [
        statements[0] == (true_count >= 3),  # Statement 1
        statements[1] == (false_count <= 5),  # Statement 2
        statements[2] == (true_count == 7),  # Statement 3
        statements[3] == (false_count == 1),  # Statement 4
        statements[4] == (statements[2] != statements[3] and (statements[2] or statements[3])),  # Statement 5
        statements[5] == is_prime(true_count),  # Statement 6
        statements[6] == is_composite(false_count)  # Statement 7
    ]
    
    # If all checks pass, this is a valid solution
    return all(checks)

# Generate all possible combinations and count valid solutions
solutions = []
for combination in product([True, False], repeat=7):
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")