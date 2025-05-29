def check_statements(true_statements):
    # Count of true statements
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Create a list of boolean values for each statement
    is_true = [i in true_statements for i in range(7)]
    
    # Check each statement
    checks = [False] * 7
    
    # Statement 1: At least 4 are true
    checks[0] = (num_true >= 4)
    
    # Statement 2: At most 3 are false
    checks[1] = (num_false <= 3)
    
    # Statement 3: Exactly 4 are true
    checks[2] = (num_true == 4)
    
    # Statement 4: Exactly 3 are false
    checks[3] = (num_false == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    checks[4] = (is_true[2] != is_true[3]) and (is_true[2] or is_true[3])
    
    # Statement 6: Number of true statements is prime
    checks[5] = num_true in [2, 3, 5, 7]
    
    # Statement 7: Number of false statements is composite
    checks[6] = num_false in [4, 6]
    
    # Check if the combination is consistent
    return all(is_true[i] == checks[i] for i in range(7))

# Generate all possible combinations and check them
solutions = []
for i in range(128):  # 2^7 possibilities
    true_statements = set()
    for j in range(7):
        if (i >> j) & 1:
            true_statements.add(j)
    if check_statements(true_statements):
        solutions.append(true_statements)

print(f"Number of solutions: {len(solutions)}")