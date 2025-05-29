def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(states):
    # Count true and false statements
    true_count = sum(states)
    false_count = 7 - true_count
    
    # Check each statement
    checks = [
        true_count >= 4,                    # Statement 1
        false_count <= 3,                   # Statement 2
        true_count == 0,                    # Statement 3
        false_count == 2,                   # Statement 4
        (states[2] ^ states[3]),           # Statement 5 (XOR)
        is_prime(true_count),              # Statement 6
        is_composite(false_count)          # Statement 7
    ]
    
    # The state of each statement must match our checks
    return all(states[i] == checks[i] for i in range(7))

# Generate and test all possible combinations
solutions = []
for i in range(2**7):
    states = [(i >> j) & 1 for j in range(7)]
    if check_statements(states):
        solutions.append(states)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"States: {sol}, True count: {sum(sol)}, False count: {7-sum(sol)}")