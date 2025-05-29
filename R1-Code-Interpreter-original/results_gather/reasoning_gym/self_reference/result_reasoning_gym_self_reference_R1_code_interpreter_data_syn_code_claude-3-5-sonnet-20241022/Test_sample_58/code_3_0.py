def verify_solution(comb):
    true_count = sum(comb)
    false_count = 7 - true_count
    
    # Verify each statement's truth value matches its claim
    
    # Statement 1: "At least 2 true"
    if comb[0] != (true_count >= 2):
        return False
        
    # Statement 2: "At most 4 false"
    if comb[1] != (false_count <= 4):
        return False
        
    # Statement 3: "Exactly 5 true"
    if comb[2] != (true_count == 5):
        return False
        
    # Statement 4: "Exactly 2 false"
    if comb[3] != (false_count == 2):
        return False
        
    # Statement 5: "Either Statement 3 or 4 is true, but not both"
    if comb[4] != (comb[2] != comb[3] and (comb[2] or comb[3])):
        return False
        
    # Statement 6: "Number of true is prime"
    if comb[5] != (true_count in {2,3,5,7}):
        return False
        
    # Statement 7: "Number of false is composite"
    if comb[6] != (false_count in {4,6}):
        return False
        
    return True

valid_solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if verify_solution(combination):
        valid_solutions.append(combination)
        print(f"Valid solution: {combination}")
        print(f"True statements: {sum(combination)}")
        print(f"False statements: {7-sum(combination)}")
        print("---")

print(f"Number of valid solutions: {len(valid_solutions)}")