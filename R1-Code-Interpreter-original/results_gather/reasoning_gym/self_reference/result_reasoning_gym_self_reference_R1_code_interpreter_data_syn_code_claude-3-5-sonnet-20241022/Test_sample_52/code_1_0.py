def analyze_solution(true_stmts):
    num_true = len(true_stmts)
    num_false = 7 - num_true
    
    print(f"\nAnalyzing solution with true statements: {sorted(list(true_stmts))}")
    print(f"Number of true statements: {num_true}")
    print(f"Number of false statements: {num_false}")
    
    # For each statement, we'll check if its presence in true_stmts matches what it claims
    
    # Statement 1: At least 3 are true
    claim1 = (num_true >= 3)
    actual1 = (0 in true_stmts)
    print(f"1. 'At least 3 true': Should be {claim1}, Is {actual1}")
    
    # Statement 2: At most 3 are false
    claim2 = (num_false <= 3)
    actual2 = (1 in true_stmts)
    print(f"2. 'At most 3 false': Should be {claim2}, Is {actual2}")
    
    # Statement 3: Exactly 4 are true
    claim3 = (num_true == 4)
    actual3 = (2 in true_stmts)
    print(f"3. 'Exactly 4 true': Should be {claim3}, Is {actual3}")
    
    # Statement 4: Exactly 3 are false
    claim4 = (num_false == 3)
    actual4 = (3 in true_stmts)
    print(f"4. 'Exactly 3 false': Should be {claim4}, Is {actual4}")
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    claim5 = (actual3 != actual4)
    actual5 = (4 in true_stmts)
    print(f"5. 'Stmt3 XOR Stmt4': Should be {claim5}, Is {actual5}")
    
    # Statement 6: Number of true statements is prime
    claim6 = (num_true in {2, 3, 5, 7})
    actual6 = (5 in true_stmts)
    print(f"6. 'Num true is prime': Should be {claim6}, Is {actual6}")
    
    # Statement 7: Number of false statements is composite
    claim7 = (num_false in {4, 6})
    actual7 = (6 in true_stmts)
    print(f"7. 'Num false is composite': Should be {claim7}, Is {actual7}")
    
    # Check if solution is valid
    valid = all(claim == actual for claim, actual in [
        (claim1, actual1), (claim2, actual2), (claim3, actual3),
        (claim4, actual4), (claim5, actual5), (claim6, actual6),
        (claim7, actual7)
    ])
    print(f"Solution is {'valid' if valid else 'invalid'}")
    return valid

solutions = [
    set(),                # Solution 1: All false
    {0, 1, 2, 3},        # Solution 2: First four true
    {6},                 # Solution 3: Only 7th true
    {0, 5, 6}           # Solution 4: 1st, 6th, and 7th true
]

valid_count = sum(analyze_solution(sol) for sol in solutions)
print(f"\nTotal number of valid solutions: {valid_count}")