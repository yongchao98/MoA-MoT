def check_statements(truth_values):
    # Count true statements
    true_count = sum(truth_values)
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
        false_count <= 1,                   # Statement 2
        true_count == 3,                    # Statement 3
        false_count == 1,                   # Statement 4
        (truth_values[2] != truth_values[3] and 
         (truth_values[2] or truth_values[3])),  # Statement 5
        is_prime(true_count),               # Statement 6
        is_composite(false_count)           # Statement 7
    ]
    
    # Check if the truth values match the statements
    return all(truth == statement for truth, statement in zip(truth_values, statements))

# Generate all possible combinations and count valid solutions
solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    truth_values = [bool((i >> j) & 1) for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {[int(x) for x in sol]}")
        print(f"Number of true statements: {sum(sol)}")