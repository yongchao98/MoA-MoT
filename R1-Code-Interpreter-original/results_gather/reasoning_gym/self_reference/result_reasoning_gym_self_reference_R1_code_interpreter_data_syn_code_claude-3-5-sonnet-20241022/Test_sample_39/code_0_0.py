def check_statements(pattern):
    # Count true and false statements
    true_count = sum(pattern)
    false_count = 7 - true_count
    
    # Check each statement
    results = [False] * 7
    
    # Statement 1: At least 4 are true
    results[0] = (true_count >= 4) == pattern[0]
    
    # Statement 2: At most 2 are false
    results[1] = (false_count <= 2) == pattern[1]
    
    # Statement 3: Exactly 1 is true
    results[2] = (true_count == 1) == pattern[2]
    
    # Statement 4: All 7 are false
    results[3] = (true_count == 0) == pattern[3]
    
    # Statement 5: Either Statement 3 or 4 is true, but not both
    stmt3_true = pattern[2]
    stmt4_true = pattern[3]
    results[4] = ((stmt3_true != stmt4_true) and (stmt3_true or stmt4_true)) == pattern[4]
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    results[5] = (true_count in primes) == pattern[5]
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    results[6] = (false_count in composites) == pattern[6]
    
    # Check if this pattern is consistent
    return all(results)

# Try all possible combinations
valid_solutions = []
for i in range(2**7):
    pattern = [bool((i >> j) & 1) for j in range(7)]
    if check_statements(pattern):
        valid_solutions.append(pattern)

print(f"Number of valid solutions: {len(valid_solutions)}")
if valid_solutions:
    print("\nValid solutions:")
    for solution in valid_solutions:
        print(f"Pattern: {[int(x) for x in solution]}, True count: {sum(solution)}")