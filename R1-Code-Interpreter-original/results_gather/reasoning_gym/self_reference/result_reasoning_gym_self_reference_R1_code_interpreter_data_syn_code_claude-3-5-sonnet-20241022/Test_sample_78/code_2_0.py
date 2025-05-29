def verify_solution(truths):
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # Verify each statement
    print(f"\nVerifying solution: {truths}")
    print(f"True count: {true_count}, False count: {false_count}")
    
    # Statement 1: "At least 7 are true"
    print(f"Statement 1: Should be {truths[0]}, 'At least 7 are true' is {true_count >= 7}")
    
    # Statement 2: "At most 5 are false"
    print(f"Statement 2: Should be {truths[1]}, 'At most 5 are false' is {false_count <= 5}")
    
    # Statement 3: "Exactly 4 are true"
    print(f"Statement 3: Should be {truths[2]}, 'Exactly 4 are true' is {true_count == 4}")
    
    # Statement 4: "Exactly 3 are false"
    print(f"Statement 4: Should be {truths[3]}, 'Exactly 3 are false' is {false_count == 3}")
    
    # Statement 5: "Either Statement 3 or 4 is true, but not both"
    s5_check = (truths[2] != truths[3]) and (truths[2] or truths[3])
    print(f"Statement 5: Should be {truths[4]}, 'Either 3 or 4 true but not both' is {s5_check}")
    
    # Statement 6: "Number of true statements is prime"
    primes = {2, 3, 5, 7}
    print(f"Statement 6: Should be {truths[5]}, 'True count is prime' is {true_count in primes}")
    
    # Statement 7: "Number of false statements is composite"
    composites = {4, 6}
    print(f"Statement 7: Should be {truths[6]}, 'False count is composite' is {false_count in composites}")

# Test all found solutions
solutions = [
    [False, False, False, False, False, False, False],
    [False, True, False, False, False, True, False],
    [False, False, False, False, False, False, True],
    [False, True, False, False, False, True, True]
]

for solution in solutions:
    verify_solution(solution)