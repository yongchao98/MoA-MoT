def verify_solution(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Create actual truth values based on the statements
    stmt1 = (num_true >= 4)
    stmt2 = (num_false <= 6)
    stmt3 = (num_true == 4)
    stmt4 = (num_false == 3)
    stmt5 = (stmt3 != stmt4)  # Either stmt3 or stmt4 is true, but not both
    stmt6 = num_true in {2,3,5,7}  # prime numbers
    stmt7 = num_false in {4,6,8,9}  # composite numbers
    
    actual_values = [stmt1, stmt2, stmt3, stmt4, stmt5, stmt6, stmt7]
    
    # Check if the proposed solution matches the actual truth values
    for i in range(7):
        if (i in true_statements) != actual_values[i]:
            return False
            
    # Additional check for logical consistency
    if stmt3 and stmt4:  # Can't both be true as they contradict
        return False
    if num_true + num_false != 7:  # Must sum to total statements
        return False
        
    return True

# Test all possible combinations with 0-7 true statements
solutions = []
for i in range(128):  # 2^7
    true_set = {j for j in range(7) if (i & (1 << j))}
    if verify_solution(true_set):
        solutions.append(true_set)

print(f"Number of valid solutions: {len(solutions)}")
for s in solutions:
    print(f"True statements: {sorted(list(s))}")