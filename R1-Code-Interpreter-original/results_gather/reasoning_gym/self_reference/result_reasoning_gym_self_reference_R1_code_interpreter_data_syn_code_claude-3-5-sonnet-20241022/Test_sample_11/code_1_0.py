def verify_solution_detailed(true_statements):
    # Convert to boolean array
    statements = [False] * 7
    for i in true_statements:
        statements[i-1] = True
    
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # For each statement, check if its truth value matches its claim
    
    # Statement 1: "At least 4 are true"
    claim1 = num_true >= 4
    valid1 = (claim1 == statements[0])
    
    # Statement 2: "At most 5 are false"
    claim2 = num_false <= 5
    valid2 = (claim2 == statements[1])
    
    # Statement 3: "Exactly 7 are true"
    claim3 = num_true == 7
    valid3 = (claim3 == statements[2])
    
    # Statement 4: "Exactly 1 is false"
    claim4 = num_false == 1
    valid4 = (claim4 == statements[3])
    
    # Statement 5: "Statement 3 XOR Statement 4"
    claim5 = (statements[2] != statements[3]) and (statements[2] or statements[3])
    valid5 = (claim5 == statements[4])
    
    # Statement 6: "Number of true statements is prime"
    claim6 = num_true in [2,3,5,7]
    valid6 = (claim6 == statements[5])
    
    # Statement 7: "Number of false statements is composite"
    claim7 = num_false in [4,6]
    valid7 = (claim7 == statements[6])
    
    print(f"\nAnalyzing solution: {true_statements}")
    print(f"True statements: {num_true}, False statements: {num_false}")
    print(f"Statement 1: {'True' if statements[0] else 'False'}, Claim valid: {valid1}")
    print(f"Statement 2: {'True' if statements[1] else 'False'}, Claim valid: {valid2}")
    print(f"Statement 3: {'True' if statements[2] else 'False'}, Claim valid: {valid3}")
    print(f"Statement 4: {'True' if statements[3] else 'False'}, Claim valid: {valid4}")
    print(f"Statement 5: {'True' if statements[4] else 'False'}, Claim valid: {valid5}")
    print(f"Statement 6: {'True' if statements[5] else 'False'}, Claim valid: {valid6}")
    print(f"Statement 7: {'True' if statements[6] else 'False'}, Claim valid: {valid7}")
    
    return all([valid1, valid2, valid3, valid4, valid5, valid6, valid7])

# Test each potential solution
solutions = [[], [2,6], [7], [2,6,7]]
valid_solutions = []

for sol in solutions:
    if verify_solution_detailed(sol):
        valid_solutions.append(sol)

print("\nFinal count of valid solutions:", len(valid_solutions))
print("Valid solutions:", valid_solutions)