def analyze_solution(true_stmts):
    num_true = len(true_stmts)
    num_false = 7 - num_false
    
    print(f"\nAnalyzing solution with true statements: {sorted(list(true_stmts))}")
    print(f"Number of true statements: {num_true}")
    print(f"Number of false statements: {7 - num_true}")
    
    # Statement 1: At least 3 are true
    stmt1_should = (num_true >= 3)
    stmt1_is = 0 in true_stmts
    print(f"Statement 1: Should be {stmt1_should}, Is {stmt1_is}")
    
    # Statement 2: At most 3 are false
    stmt2_should = ((7 - num_true) <= 3)
    stmt2_is = 1 in true_stmts
    print(f"Statement 2: Should be {stmt2_should}, Is {stmt2_is}")
    
    # Statement 3: Exactly 4 are true
    stmt3_should = (num_true == 4)
    stmt3_is = 2 in true_stmts
    print(f"Statement 3: Should be {stmt3_should}, Is {stmt3_is}")
    
    # Statement 4: Exactly 3 are false
    stmt4_should = ((7 - num_true) == 3)
    stmt4_is = 3 in true_stmts
    print(f"Statement 4: Should be {stmt4_should}, Is {stmt4_is}")
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    stmt5_should = ((2 in true_stmts) != (3 in true_stmts))
    stmt5_is = 4 in true_stmts
    print(f"Statement 5: Should be {stmt5_should}, Is {stmt5_is}")
    
    # Statement 6: Number of true statements is prime
    stmt6_should = (num_true in {2, 3, 5, 7})
    stmt6_is = 5 in true_stmts
    print(f"Statement 6: Should be {stmt6_should}, Is {stmt6_is}")
    
    # Statement 7: Number of false statements is composite
    stmt7_should = ((7 - num_true) in {4, 6})
    stmt7_is = 6 in true_stmts
    print(f"Statement 7: Should be {stmt7_should}, Is {stmt7_is}")

solutions = [
    set(),                # Solution 1: All false
    {0, 1, 2, 3},        # Solution 2: First four true
    {6},                 # Solution 3: Only 7th true
    {0, 5, 6}           # Solution 4: 1st, 6th, and 7th true
]

for sol in solutions:
    analyze_solution(sol)