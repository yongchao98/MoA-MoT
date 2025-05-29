def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def verify_solution_detailed(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    print(f"\nDetailed verification for solution: {sorted(list(true_statements))}")
    print(f"Number of true statements: {num_true}")
    print(f"Number of false statements: {num_false}")
    
    # Verify each statement's actual truth value
    stmt1_actual = num_true >= 3
    stmt2_actual = num_false <= 3
    stmt3_actual = num_true == 4
    stmt4_actual = num_false == 3
    stmt5_actual = (2 in true_statements) != (3 in true_statements)
    stmt6_actual = is_prime(num_true)
    stmt7_actual = is_composite(num_false)
    
    actuals = [stmt1_actual, stmt2_actual, stmt3_actual, stmt4_actual, 
               stmt5_actual, stmt6_actual, stmt7_actual]
    claimed = [i in true_statements for i in range(7)]
    
    all_valid = True
    for i, (actual, claim) in enumerate(zip(actuals, claimed)):
        if actual != claim:
            all_valid = False
        print(f"Statement {i+1}:")
        print(f"  Actual value: {actual}")
        print(f"  Claimed value: {claim}")
        print(f"  Match: {actual == claim}")
    
    print(f"\nSolution is {'valid' if all_valid else 'invalid'}")
    return all_valid

# Verify the two solutions we've seen are valid
solutions = [set(), {0, 1, 2, 3}]
valid_count = sum(1 for s in solutions if verify_solution_detailed(s))
print(f"\nTotal confirmed valid solutions: {valid_count}")