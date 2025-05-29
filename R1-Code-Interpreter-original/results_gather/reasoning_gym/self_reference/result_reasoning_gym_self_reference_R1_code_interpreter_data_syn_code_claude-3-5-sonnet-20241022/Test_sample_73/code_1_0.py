def verify_solution(solution_set):
    # Convert to list of booleans
    statements = [i in solution_set for i in range(7)]
    num_true = sum(statements)
    num_false = 7 - num_true
    
    # Verify each statement
    # Statement 1: At least 5 are true
    s1 = (num_true >= 5) == statements[0]
    
    # Statement 2: At most 6 are false
    s2 = (num_false <= 6) == statements[1]
    
    # Statement 3: Exactly 5 are true
    s3 = (num_true == 5) == statements[2]
    
    # Statement 4: Exactly 2 are false
    s4 = (num_false == 2) == statements[3]
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    s5 = ((statements[2] != statements[3]) and (statements[2] or statements[3])) == statements[4]
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    s6 = (num_true in primes) == statements[5]
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}
    s7 = (num_false in composites) == statements[6]
    
    return all([s1, s2, s3, s4, s5, s6, s7])

solutions = [
    set(),
    {1, 5},
    {0, 1, 2, 3, 5},
    {1, 5, 6}
]

print("Verification of solutions:")
for i, sol in enumerate(solutions):
    print(f"Solution {i+1} ({sorted(list(sol))}): {verify_solution(sol)}")