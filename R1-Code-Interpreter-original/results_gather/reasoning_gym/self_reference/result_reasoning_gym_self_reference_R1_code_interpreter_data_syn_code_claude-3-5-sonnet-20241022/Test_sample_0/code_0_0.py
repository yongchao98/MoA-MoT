def check_statements(true_statements):
    # Count of true statements
    num_true = len(true_statements)
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
    results = [False] * 7
    
    # Statement 1: At least 3 are true
    results[0] = (num_true >= 3)
    
    # Statement 2: At most 1 is false
    results[1] = (num_false <= 1)
    
    # Statement 3: Exactly 3 are true
    results[2] = (num_true == 3)
    
    # Statement 4: Exactly 5 are false
    results[3] = (num_false == 5)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = ((3 in true_statements) != (4 in true_statements))
    
    # Statement 6: Number of true statements is prime
    results[5] = is_prime(num_true)
    
    # Statement 7: Number of false statements is composite
    results[6] = is_composite(num_false)
    
    # Check if the results match the given true_statements
    return set(true_statements) == {i for i in range(7) if results[i]}

# Generate all possible combinations
valid_solutions = []
for i in range(128):  # 2^7 = 128 combinations
    true_statements = {j for j in range(7) if (i & (1 << j))}
    if check_statements(true_statements):
        valid_solutions.append(true_statements)

print(len(valid_solutions))