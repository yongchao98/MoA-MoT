def check_statements(truths):
    # Count number of true statements
    true_count = sum(truths)
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
    checks = [
        true_count >= 4,                    # Statement 1
        false_count <= 2,                   # Statement 2
        true_count == 5,                    # Statement 3
        false_count == 6,                   # Statement 4
        (truths[2] != truths[3] and (truths[2] or truths[3])),  # Statement 5
        is_prime(true_count),               # Statement 6
        is_composite(false_count)           # Statement 7
    ]
    
    # Check if the truth values match the actual evaluation of statements
    return all(truths[i] == checks[i] for i in range(7))

# Generate all possible combinations and count valid solutions
valid_solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        valid_solutions.append(truths)

print(f"Number of valid solutions: {len(valid_solutions)}")
if valid_solutions:
    print("\nValid solutions:")
    for solution in valid_solutions:
        print(f"Solution: {solution}")
        print(f"True statements: {sum(solution)}")
        print(f"False statements: {7-sum(solution)}")
        print()