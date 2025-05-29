def verify_solution(solution):
    true_count = sum(solution)
    false_count = 7 - true_count
    
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0: return False
        return True
    
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    # Check each statement's claim against the actual truth value
    checks = [
        (true_count >= 3) == solution[0],
        (false_count <= 3) == solution[1],
        (true_count == 5) == solution[2],
        (false_count == 1) == solution[3],
        (solution[2] != solution[3]) == solution[4],
        is_prime(true_count) == solution[5],
        is_composite(false_count) == solution[6]
    ]
    
    return all(checks)

solutions = [
    [False, False, False, False, False, False, False],
    [True, True, True, False, True, True, False],
    [False, False, False, False, False, False, True],
    [True, False, False, False, False, True, True]
]

print("Verifying each solution:")
for i, sol in enumerate(solutions, 1):
    print(f"Solution {i} is {'valid' if verify_solution(sol) else 'invalid'}")