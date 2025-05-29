def verify_solution(solution):
    true_count = sum(solution)
    false_count = 7 - true_count
    
    # Verify each statement
    checks = [
        # Statement 1: At least 1 statement is true
        (true_count >= 1) == solution[0],
        
        # Statement 2: At most 3 are false
        (false_count <= 3) == solution[1],
        
        # Statement 3: Exactly 4 are true
        (true_count == 4) == solution[2],
        
        # Statement 4: Exactly 3 are false
        (false_count == 3) == solution[3],
        
        # Statement 5: Statement 3 XOR Statement 4 is true
        (solution[2] != solution[3] and (solution[2] or solution[3])) == solution[4],
        
        # Statement 6: Number of true statements is prime
        (true_count in [2, 3, 5, 7]) == solution[5],
        
        # Statement 7: Number of false statements is composite
        (false_count in [4, 6]) == solution[6]
    ]
    
    return all(checks)

solutions = [
    [0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0, 1, 1]
]

print("Verifying each solution:")
for i, sol in enumerate(solutions, 1):
    print(f"Solution {i} is {'valid' if verify_solution(sol) else 'invalid'}")