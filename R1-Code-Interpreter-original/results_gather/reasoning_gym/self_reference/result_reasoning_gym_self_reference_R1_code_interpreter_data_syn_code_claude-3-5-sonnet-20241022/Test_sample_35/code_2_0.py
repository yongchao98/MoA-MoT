def validate_solution(comb):
    true_count = sum(comb)
    false_count = 7 - true_count
    
    # For each statement, check if its truth value matches what it claims
    
    # Statement 1: "At least 2 true"
    if comb[0] == True and true_count < 2: return False
    if comb[0] == False and true_count >= 2: return False
    
    # Statement 2: "At most 3 false"
    if comb[1] == True and false_count > 3: return False
    if comb[1] == False and false_count <= 3: return False
    
    # Statement 3: "Exactly 4 true"
    if comb[2] == True and true_count != 4: return False
    if comb[2] == False and true_count == 4: return False
    
    # Statement 4: "Exactly 3 false"
    if comb[3] == True and false_count != 3: return False
    if comb[3] == False and false_count == 3: return False
    
    # Statement 5: "Statement 3 XOR 4 is true"
    stmt3_true = (true_count == 4)
    stmt4_true = (false_count == 3)
    if comb[4] == True and not (stmt3_true != stmt4_true and (stmt3_true or stmt4_true)): return False
    if comb[4] == False and (stmt3_true != stmt4_true and (stmt3_true or stmt4_true)): return False
    
    # Statement 6: "True count is prime"
    primes = {2, 3, 5, 7}
    if comb[5] == True and true_count not in primes: return False
    if comb[5] == False and true_count in primes: return False
    
    # Statement 7: "False count is composite"
    composites = {4, 6}
    if comb[6] == True and false_count not in composites: return False
    if comb[6] == False and false_count in composites: return False
    
    return True

# Test all possible combinations
valid_solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if validate_solution(combination):
        valid_solutions.append(combination)
        print(f"Valid solution: {combination}, True count: {sum(combination)}, False count: {7-sum(combination)}")

print(f"\nTotal valid solutions: {len(valid_solutions)}")