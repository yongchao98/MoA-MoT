def check_statements(combination):
    # Count true statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Helper function to check if a number is prime
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    
    # Helper function to check if a number is composite
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    # Check each statement
    statements = [
        true_count >= 6,                    # Statement 1
        false_count <= 1,                   # Statement 2
        true_count == 6,                    # Statement 3
        false_count == 5,                   # Statement 4
        (combination[2] != combination[3]),  # Statement 5
        is_prime(true_count),               # Statement 6
        is_composite(false_count)           # Statement 7
    ]
    
    # Check if the combination is consistent
    return sum([s == c for s, c in zip(statements, combination)]) == true_count

# Generate and test all possible combinations
solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")
if len(solutions) > 0:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Solution: {sol} (True count: {sum(sol)})")