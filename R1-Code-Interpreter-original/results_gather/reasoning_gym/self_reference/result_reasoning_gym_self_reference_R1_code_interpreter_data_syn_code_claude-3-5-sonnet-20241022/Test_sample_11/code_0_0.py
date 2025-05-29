def check_statements(true_statements):
    # Count number of true and false statements
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Initialize all statements as False
    statements = [False] * 7
    
    # Set the specified statements as True
    for i in true_statements:
        statements[i] = True
    
    # Check each statement's validity
    valid = [
        # Statement 1: At least 4 are true
        (num_true >= 4) == statements[0],
        
        # Statement 2: At most 5 are false
        (num_false <= 5) == statements[1],
        
        # Statement 3: Exactly 7 are true
        (num_true == 7) == statements[2],
        
        # Statement 4: Exactly 1 is false
        (num_false == 1) == statements[3],
        
        # Statement 5: Statement 3 or 4 is true (XOR)
        ((statements[2] != statements[3]) and (statements[2] or statements[3])) == statements[4],
        
        # Statement 6: Number of true statements is prime
        (num_true in [2,3,5,7]) == statements[5],
        
        # Statement 7: Number of false statements is composite
        (num_false in [4,6]) == statements[6]
    ]
    
    # Return True if all statements are consistent
    return all(valid)

# Generate all possible combinations and check them
solutions = []
for i in range(128):  # 2^7 possibilities
    true_statements = []
    for j in range(7):
        if (i >> j) & 1:
            true_statements.append(j)
    
    if check_statements(true_statements):
        solutions.append(true_statements)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("Solutions found:")
    for sol in solutions:
        print(f"True statements: {[x+1 for x in sol]}")