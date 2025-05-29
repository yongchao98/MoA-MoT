def is_valid_solution(true_statements):
    # Convert to boolean array
    statements = [False] * 7
    for i in true_statements:
        statements[i-1] = True
    
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Check each statement's claim against its truth value
    
    # Statement 1: "At least 4 are true"
    if (num_true >= 4) != statements[0]:
        return False
    
    # Statement 2: "At most 5 are false"
    if (num_false <= 5) != statements[1]:
        return False
    
    # Statement 3: "Exactly 7 are true"
    if (num_true == 7) != statements[2]:
        return False
    
    # Statement 4: "Exactly 1 is false"
    if (num_false == 1) != statements[3]:
        return False
    
    # Statement 5: "Statement 3 XOR Statement 4"
    if ((statements[2] != statements[3]) and (statements[2] or statements[3])) != statements[4]:
        return False
    
    # Statement 6: "Number of true statements is prime"
    if (num_true in [2,3,5,7]) != statements[5]:
        return False
    
    # Statement 7: "Number of false statements is composite"
    if (num_false in [4,6]) != statements[6]:
        return False
    
    return True

# Test all claimed solutions
solutions = [[], [2,6], [7], [2,6,7]]
valid_solutions = []

print("Checking each solution:")
for sol in solutions:
    if is_valid_solution(sol):
        valid_solutions.append(sol)
        print(f"Solution {sol} is valid")
    else:
        print(f"Solution {sol} is invalid")

print(f"\nTotal valid solutions: {len(valid_solutions)}")
print("Valid solutions:", valid_solutions)