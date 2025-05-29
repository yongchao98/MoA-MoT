def check_statements(truths):
    # Count true and false statements
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # Check each statement
    checks = [
        # Statement 1: At least 7 statements are true
        (truths[0] == (true_count >= 7)),
        
        # Statement 2: At most 1 statement is false
        (truths[1] == (false_count <= 1)),
        
        # Statement 3: Exactly 2 statements are true
        (truths[2] == (true_count == 2)),
        
        # Statement 4: Exactly 4 statements are false
        (truths[3] == (false_count == 4)),
        
        # Statement 5: Either Statement 3 or Statement 4 is true, but not both
        (truths[4] == (truths[2] != truths[3] and (truths[2] or truths[3]))),
        
        # Statement 6: Number of true statements is prime
        (truths[5] == (true_count in [2, 3, 5, 7])),
        
        # Statement 7: Number of false statements is composite
        (truths[6] == (false_count in [4, 6]))
    ]
    
    # All checks must be True for a valid solution
    return all(checks)

# Try all possible combinations
solutions = []
for i in range(2**7):  # 2^7 = 128 possibilities
    # Convert number to binary representation of 7 bits
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truths):
        solutions.append(truths)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {sol}")
        print(f"Number of true statements: {sum(sol)}")
        print(f"Number of false statements: {7-sum(sol)}")
        print()