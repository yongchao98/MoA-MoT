def check_statements(combination):
    # Count true and false statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement's consistency
    # Statement 1: At least 2 are true
    stmt1_valid = (combination[0] == (true_count >= 2))
    
    # Statement 2: At most 4 are false
    stmt2_valid = (combination[1] == (false_count <= 4))
    
    # Statement 3: Exactly 5 are true
    stmt3_valid = (combination[2] == (true_count == 5))
    
    # Statement 4: Exactly 7 are false
    stmt4_valid = (combination[3] == (false_count == 7))
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    stmt5_valid = (combination[4] == (combination[2] != combination[3] and (combination[2] or combination[3])))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    stmt6_valid = (combination[5] == (true_count in primes))
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    stmt7_valid = (combination[6] == (false_count in composites))
    
    return all([stmt1_valid, stmt2_valid, stmt3_valid, stmt4_valid, stmt5_valid, stmt6_valid, stmt7_valid])

# Generate all possible combinations and check each
solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    combination = [bool(i & (1 << j)) for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {[int(x) for x in sol]}")
        print(f"True statements: {sum(sol)}")
        print(f"False statements: {7-sum(sol)}")
        print()