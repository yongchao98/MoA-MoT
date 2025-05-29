def verify_final(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # What each statement actually means
    actual_values = [
        num_true >= 3,                    # Statement 1
        num_false <= 3,                   # Statement 2
        num_true == 4,                    # Statement 3
        num_false == 3,                   # Statement 4
        (2 in true_statements) != (3 in true_statements),  # Statement 5
        num_true in {2,3,5,7},           # Statement 6
        num_false in {4,6}               # Statement 7
    ]
    
    # What each statement is claimed to be
    claimed_values = [i in true_statements for i in range(7)]
    
    print(f"Solution: {sorted(list(true_statements))}")
    print(f"True statements: {num_true}, False statements: {num_false}")
    
    for i, (actual, claimed) in enumerate(zip(actual_values, claimed_values)):
        if actual != claimed:
            print(f"Mismatch in statement {i+1}: actual={actual}, claimed={claimed}")
            return False
    return True

# Test solution 1: all false
print("\nTesting solution 1:")
sol1_valid = verify_final(set())

# Test solution 2: statements 1,2,3,4 true
print("\nTesting solution 2:")
sol2_valid = verify_final({0,1,2,3})

print(f"\nValid solutions found: {sol1_valid + sol2_valid}")