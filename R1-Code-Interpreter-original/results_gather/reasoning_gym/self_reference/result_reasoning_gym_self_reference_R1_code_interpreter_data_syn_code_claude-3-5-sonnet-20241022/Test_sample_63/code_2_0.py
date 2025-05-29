def verify_solution(combination):
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Verify each statement's claim against the actual state
    
    # Statement 1: "At least 2 are true"
    if combination[0] == 1 and true_count < 2:
        return False
    if combination[0] == 0 and true_count >= 2:
        return False
    
    # Statement 2: "At most 2 are false"
    if combination[1] == 1 and false_count > 2:
        return False
    if combination[1] == 0 and false_count <= 2:
        return False
    
    # Statement 3: "Exactly 4 are true"
    if combination[2] == 1 and true_count != 4:
        return False
    if combination[2] == 0 and true_count == 4:
        return False
    
    # Statement 4: "Exactly 2 are false"
    if combination[3] == 1 and false_count != 2:
        return False
    if combination[3] == 0 and false_count == 2:
        return False
    
    # Statement 5: "Either Statement 3 or Statement 4 is true, but not both"
    stmt5_condition = (combination[2] != combination[3] and (combination[2] or combination[3]))
    if combination[4] == 1 and not stmt5_condition:
        return False
    if combination[4] == 0 and stmt5_condition:
        return False
    
    # Statement 6: "Number of true statements is prime"
    primes = {2, 3, 5, 7}
    if combination[5] == 1 and true_count not in primes:
        return False
    if combination[5] == 0 and true_count in primes:
        return False
    
    # Statement 7: "Number of false statements is composite"
    composites = {4, 6}
    if combination[6] == 1 and false_count not in composites:
        return False
    if combination[6] == 0 and false_count in composites:
        return False
    
    return True

# Test all solutions
solutions = [
    [0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 1, 0],
    [1, 1, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0, 1],
    [1, 0, 0, 0, 0, 1, 1]
]

valid_solutions = 0
for i, sol in enumerate(solutions):
    if verify_solution(sol):
        valid_solutions += 1
        print(f"Solution {i+1} is valid")
    else:
        print(f"Solution {i+1} is invalid")

print(f"\nTotal valid solutions: {valid_solutions}")