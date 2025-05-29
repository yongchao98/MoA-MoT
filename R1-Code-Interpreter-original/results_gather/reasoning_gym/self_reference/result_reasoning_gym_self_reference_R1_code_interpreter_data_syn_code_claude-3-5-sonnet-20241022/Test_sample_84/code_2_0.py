def verify_single_solution(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Evaluate what each statement claims
    stmt1 = (num_true >= 3)                    # At least 3 true
    stmt2 = (num_false <= 3)                   # At most 3 false
    stmt3 = (num_true == 4)                    # Exactly 4 true
    stmt4 = (num_false == 3)                   # Exactly 3 false
    stmt5 = ((2 in true_statements) != (3 in true_statements))  # Statement 3 XOR 4
    stmt6 = (num_true in {2,3,5,7})           # Num true is prime
    stmt7 = (num_false in {4,6})              # Num false is composite
    
    actual_values = [stmt1, stmt2, stmt3, stmt4, stmt5, stmt6, stmt7]
    claimed_values = [i in true_statements for i in range(7)]
    
    matches = [actual == claimed for actual, claimed in zip(actual_values, claimed_values)]
    
    print(f"Solution with true statements: {sorted(list(true_statements))}")
    print(f"True statements: {num_true}, False statements: {num_false}")
    print("Statement evaluation:")
    for i, (actual, claimed, match) in enumerate(zip(actual_values, claimed_values, matches)):
        print(f"Statement {i+1}: Actual={actual}, Claimed={claimed}, Match={match}")
    print(f"Valid: {all(matches)}\n")
    
    return all(matches)

# Test all solutions from previous output
solutions = [
    set(),               # no true statements
    {0, 1, 2, 3},       # statements 1,2,3,4 are true
    {6},                # only statement 7 is true
    {0, 5, 6}          # statements 1,6,7 are true
]

valid_solutions = []
for sol in solutions:
    if verify_single_solution(sol):
        valid_solutions.append(sol)

print(f"Number of valid solutions: {len(valid_solutions)}")