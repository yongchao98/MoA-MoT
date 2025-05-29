def check_consistency(true_count):
    # Given the number of true statements, check if it's consistent
    false_count = 7 - true_count
    
    # Initialize all statements as false
    statements = [False] * 7
    
    # For each statement, determine if it should be true based on true_count
    
    # Statement 1: "At least 7 of these 7 statements are true"
    stmt1 = (true_count >= 7)
    
    # Statement 2: "At most 7 of these 7 statements are false"
    stmt2 = (false_count <= 7)  # Always true
    
    # Statement 3: "Exactly 5 of these 7 statements are true"
    stmt3 = (true_count == 5)
    
    # Statement 4: "Exactly 0 of these 7 statements are false"
    stmt4 = (false_count == 0)
    
    # Statement 5: "Either Statement 3 or Statement 4 is true, but not both"
    stmt5 = (stmt3 != stmt4) and (stmt3 or stmt4)
    
    # Statement 6: "The number of true statements is a prime number"
    primes = {2, 3, 5, 7}
    stmt6 = (true_count in primes)
    
    # Statement 7: "The number of false statements is a composite number"
    composites = {4, 6}
    stmt7 = (false_count in composites)
    
    # Collect all statements
    statements = [stmt1, stmt2, stmt3, stmt4, stmt5, stmt6, stmt7]
    
    # Count true statements in our evaluation
    actual_true_count = sum(statements)
    
    # Solution is valid if the number of true statements matches our assumption
    return actual_true_count == true_count, statements

# Check all possible numbers of true statements (0 to 7)
solutions = []
for true_count in range(8):
    is_valid, statements = check_consistency(true_count)
    if is_valid:
        solutions.append((true_count, statements))

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for true_count, statements in solutions:
        print(f"Number of true statements: {true_count}")
        print(f"Truth values: {statements}")
        print(f"Number of false statements: {7-true_count}")
        print()