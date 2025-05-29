def verify_solution(true_stmts):
    num_true = len(true_stmts)
    num_false = 7 - num_true
    
    # For each statement, check if it's marked as true and if it should be true
    checks = [
        (0, num_true >= 3),                    # Statement 1
        (1, num_false <= 3),                   # Statement 2
        (2, num_true == 4),                    # Statement 3
        (3, num_false == 3),                   # Statement 4
        (4, (2 in true_stmts) != (3 in true_stmts)),  # Statement 5
        (5, num_true in {2, 3, 5, 7}),        # Statement 6
        (6, num_false in {4, 6})              # Statement 7
    ]
    
    # Check if each statement's truth value matches what it should be
    for idx, should_be_true in checks:
        is_true = idx in true_stmts
        if is_true != should_be_true:
            return False
    return True

solutions = [
    set(),
    {0, 1, 2, 3},
    {6},
    {0, 5, 6}
]

valid_solutions = []
for sol in solutions:
    if verify_solution(sol):
        valid_solutions.append(sol)

print(f"Number of valid solutions: {len(valid_solutions)}")
print("\nValid solutions:")
for sol in valid_solutions:
    print(f"True statements: {sorted(list(sol))}")