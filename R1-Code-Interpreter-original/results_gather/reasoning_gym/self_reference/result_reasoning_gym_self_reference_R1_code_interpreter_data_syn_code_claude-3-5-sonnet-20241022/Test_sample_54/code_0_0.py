def check_statements(truths):
    # Count number of true statements
    num_true = sum(truths)
    num_false = 7 - num_true
    
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
    results = [
        num_true >= 6,                    # Statement 1
        num_false <= 2,                   # Statement 2
        num_true == 6,                    # Statement 3
        num_false == 0,                   # Statement 4
        (truths[2] != truths[3] and      # Statement 5
         (truths[2] or truths[3])),
        is_prime(num_true),              # Statement 6
        is_composite(num_false)          # Statement 7
    ]
    
    # Check if the assignment of truth values matches the results
    return all(truths[i] == results[i] for i in range(7))

# Generate all possible combinations and check them
valid_solutions = []
for i in range(2**7):
    # Convert number to binary representation of true/false values
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        valid_solutions.append(truths)

print(f"Number of valid solutions: {len(valid_solutions)}")
if valid_solutions:
    print("\nValid solutions:")
    for solution in valid_solutions:
        print(solution)