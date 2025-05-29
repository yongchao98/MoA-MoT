def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def verify_single_solution(truths):
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # Verify each statement
    checks = [
        # Statement 1: At least 7 are true
        truths[0] == (true_count >= 7),
        
        # Statement 2: At most 5 are false
        truths[1] == (false_count <= 5),
        
        # Statement 3: Exactly 4 are true
        truths[2] == (true_count == 4),
        
        # Statement 4: Exactly 3 are false
        truths[3] == (false_count == 3),
        
        # Statement 5: Either Statement 3 or 4 is true, but not both
        truths[4] == ((truths[2] != truths[3]) and (truths[2] or truths[3])),
        
        # Statement 6: True count is prime
        truths[5] == is_prime(true_count),
        
        # Statement 7: False count is composite
        truths[6] == (not is_prime(false_count) and false_count > 1)
    ]
    
    return all(checks)

# Test all claimed solutions
solutions = [
    [False, False, False, False, False, False, False],
    [False, True, False, False, False, True, False],
    [False, False, False, False, False, False, True],
    [False, True, False, False, False, True, True]
]

valid_solutions = [s for s in solutions if verify_single_solution(s)]
print(f"Number of valid solutions: {len(valid_solutions)}")