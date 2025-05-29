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
    statements = [
        true_count >= 5,                    # Statement 1
        false_count <= 3,                   # Statement 2
        true_count == 5,                    # Statement 3
        false_count == 5,                   # Statement 4
        (truths[2] != truths[3]),          # Statement 5
        is_prime(true_count),              # Statement 6
        is_composite(false_count)          # Statement 7
    ]
    
    # Check if the actual truth values match the computed statements
    return all(truths[i] == statements[i] for i in range(7))

# Generate all possible combinations and check each
solutions = []
for i in range(2**7):
    # Convert number to binary representation of true/false values
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        solutions.append(truths)
        print(f"Solution found: {truths}")
        print(f"Number of true statements: {sum(truths)}")
        print(f"Number of false statements: {7-sum(truths)}")
        print("---")

print(f"Total number of solutions: {len(solutions)}")