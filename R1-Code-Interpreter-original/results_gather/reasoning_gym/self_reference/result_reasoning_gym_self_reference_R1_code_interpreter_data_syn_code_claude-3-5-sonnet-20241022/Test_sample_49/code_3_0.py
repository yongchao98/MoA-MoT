def verify_solution(combination):
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # For each statement, check if its truth value matches what it claims
    
    # Statement 1: "At least 7 are true"
    if combination[0] == True and true_count < 7:
        return False
    if combination[0] == False and true_count >= 7:
        return False
    
    # Statement 2: "At most 2 are false"
    if combination[1] == True and false_count > 2:
        return False
    if combination[1] == False and false_count <= 2:
        return False
    
    # Statement 3: "Exactly 5 are true"
    if combination[2] == True and true_count != 5:
        return False
    if combination[2] == False and true_count == 5:
        return False
    
    # Statement 4: "Exactly 2 are false"
    if combination[3] == True and false_count != 2:
        return False
    if combination[3] == False and false_count == 2:
        return False
    
    # Statement 5: "Either Statement 3 or Statement 4 is true, but not both"
    stmt5_should_be = (combination[2] != combination[3]) and (combination[2] or combination[3])
    if combination[4] != stmt5_should_be:
        return False
    
    # Statement 6: "Number of true statements is prime"
    primes = {2, 3, 5, 7}
    if combination[5] == True and true_count not in primes:
        return False
    if combination[5] == False and true_count in primes:
        return False
    
    # Statement 7: "Number of false statements is composite"
    composites = {4, 6}  # Only relevant composites for this problem
    if combination[6] == True and false_count not in composites:
        return False
    if combination[6] == False and false_count in composites:
        return False
    
    return True

# Test all possible combinations
valid_solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if verify_solution(combination):
        valid_solutions.append(combination)

print(f"Number of valid solutions: {len(valid_solutions)}")
for sol in valid_solutions:
    print(f"\nSolution: {sol}")
    print(f"True statements: {sum(sol)}")
    print(f"False statements: {7-sum(sol)}")
    print("Statement values:")
    for i, val in enumerate(sol, 1):
        print(f"Statement {i}: {'True' if val else 'False'}")