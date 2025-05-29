def check_statements(values):
    # Count true and false statements
    true_count = sum(values)
    false_count = 7 - true_count
    
    # Check each statement
    s1 = (true_count >= 1) == values[0]  # At least 1 true
    s2 = (false_count <= 2) == values[1]  # At most 2 false
    s3 = (true_count == 0) == values[2]  # Exactly 0 true
    s4 = (false_count == 0) == values[3]  # Exactly 0 false
    s5 = ((values[2] != values[3]) and (not (values[2] and values[3]))) == values[4]  # XOR of s3 and s4
    
    # Check if true_count is prime
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0: return False
        return True
    
    # Check if false_count is composite
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    s6 = is_prime(true_count) == values[5]  # True count is prime
    s7 = is_composite(false_count) == values[6]  # False count is composite
    
    return all([s1, s2, s3, s4, s5, s6, s7])

# Generate all possible combinations and check each
solutions = []
for i in range(2**7):
    values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(values):
        solutions.append(values)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Values: {sol}, True count: {sum(sol)}, False count: {7-sum(sol)}")