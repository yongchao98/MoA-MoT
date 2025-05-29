def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(truth_values):
    # Count true statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 7 statements are true
    results[0] = (true_count >= 7)
    
    # Statement 2: At most 2 statements are false
    results[1] = (false_count <= 2)
    
    # Statement 3: Exactly 7 statements are true
    results[2] = (true_count == 7)
    
    # Statement 4: Exactly 5 statements are false
    results[3] = (false_count == 5)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    results[4] = (results[2] != results[3]) and (results[2] or results[3])
    
    # Statement 6: Number of true statements is prime
    results[5] = is_prime(true_count)
    
    # Statement 7: Number of false statements is composite
    results[6] = is_composite(false_count)
    
    # Check if this combination is valid
    return list(results) == list(truth_values)

# Generate and test all possible combinations
valid_solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truth_values):
        valid_solutions.append(truth_values)

print(f"Number of valid solutions: {len(valid_solutions)}")
if valid_solutions:
    print("\nValid solutions:")
    for solution in valid_solutions:
        print(solution)