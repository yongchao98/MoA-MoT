def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_consistency(truths):
    # Get initial counts
    true_count = sum(truths)
    false_count = 7 - true_count
    
    # Check each statement
    checks = [
        true_count >= 4,              # Statement 1
        false_count <= 1,             # Statement 2
        true_count == 4,              # Statement 3
        false_count == 3,             # Statement 4
        (truths[2] ^ truths[3]),      # Statement 5 (XOR of statements 3 and 4)
        is_prime(true_count),         # Statement 6
        is_composite(false_count)     # Statement 7
    ]
    
    # If the checks match the assumed truths, it's consistent
    return all(check == truth for check, truth in zip(checks, truths))

# Try all possibilities
solutions = []
for i in range(2**7):
    truths = [(i >> j) & 1 == 1 for j in range(7)]
    if check_consistency(truths):
        solutions.append(truths)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    for i, sol in enumerate(solutions, 1):
        print(f"\nSolution {i}:")
        print(f"Truth values: {sol}")
        print(f"True statements: {sum(sol)}")
        print(f"False statements: {7-sum(sol)}")