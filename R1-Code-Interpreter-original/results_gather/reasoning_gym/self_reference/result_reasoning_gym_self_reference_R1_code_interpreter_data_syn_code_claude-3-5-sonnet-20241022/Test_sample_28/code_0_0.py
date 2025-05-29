def check_statements(true_statements):
    # Count of true statements
    true_count = len(true_statements)
    false_count = 7 - true_count
    
    # Create a list of boolean values for each statement
    is_true = [i+1 in true_statements for i in range(7)]
    
    # Check each statement
    checks = [
        # Statement 1: At least 5 statements are true
        (true_count >= 5) == is_true[0],
        
        # Statement 2: At most 2 statements are false
        (false_count <= 2) == is_true[1],
        
        # Statement 3: Exactly 3 statements are true
        (true_count == 3) == is_true[2],
        
        # Statement 4: Exactly 1 statement is false
        (false_count == 1) == is_true[3],
        
        # Statement 5: Either Statement 3 or Statement 4 is true, but not both
        ((is_true[2] != is_true[3]) and (is_true[2] or is_true[3])) == is_true[4],
        
        # Statement 6: Number of true statements is prime
        (true_count in [2,3,5,7]) == is_true[5],
        
        # Statement 7: Number of false statements is composite
        (false_count in [4,6]) == is_true[6]
    ]
    
    # All checks must be True for a valid solution
    return all(checks)

# Generate all possible combinations
solutions = []
for i in range(128):  # 2^7 possibilities
    true_statements = set()
    for j in range(7):
        if (i >> j) & 1:
            true_statements.add(j+1)
    if check_statements(true_statements):
        solutions.append(true_statements)

print(f"Number of solutions: {len(solutions)}")
print("\nSolutions:")
for sol in solutions:
    print(f"True statements: {sorted(sol)}")