def check_statements(truth_values):
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # For each truth assignment, calculate what each statement SHOULD be
    should_be = [False] * 7
    
    # Statement 1: At least 7 of these 7 statements are true
    should_be[0] = (true_count >= 7)
    
    # Statement 2: At most 4 of these 7 statements are false
    should_be[1] = (false_count <= 4)
    
    # Statement 3: Exactly 0 of these 7 statements are true
    should_be[2] = (true_count == 0)
    
    # Statement 4: Exactly 3 of these 7 statements are false
    should_be[3] = (false_count == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    should_be[4] = bool(truth_values[2]) != bool(truth_values[3])
    
    # Statement 6: The number of true statements is a prime number
    primes = {2, 3, 5, 7}
    should_be[5] = (true_count in primes)
    
    # Statement 7: The number of false statements is a composite number
    composites = {4, 6}
    should_be[6] = (false_count in composites)
    
    # Check if the actual truth values match what they should be
    return should_be == list(truth_values)

# Test all possible combinations
solutions = []
for i in range(2**7):
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {sol}")
        print(f"Number of true statements: {sum(sol)}")
        print(f"Number of false statements: {7-sum(sol)}")
        print()