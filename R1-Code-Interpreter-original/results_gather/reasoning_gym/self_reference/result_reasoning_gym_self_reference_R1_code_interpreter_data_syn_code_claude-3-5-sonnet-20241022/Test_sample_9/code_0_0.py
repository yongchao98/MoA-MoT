def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_combination(comb):
    # Count true statements
    true_count = sum(comb)
    false_count = 7 - true_count
    
    # Check each statement
    checks = [False] * 7
    
    # Statement 1: At least 3 are true
    checks[0] = (true_count >= 3)
    
    # Statement 2: At most 7 are false
    checks[1] = (false_count <= 7)
    
    # Statement 3: Exactly 4 are true
    checks[2] = (true_count == 4)
    
    # Statement 4: Exactly 3 are false
    checks[3] = (false_count == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    checks[4] = (comb[2] != comb[3] and (comb[2] or comb[3]))
    
    # Statement 6: Number of true statements is prime
    checks[5] = is_prime(true_count)
    
    # Statement 7: Number of false statements is composite
    checks[6] = is_composite(false_count)
    
    # Check if the combination matches the actual true/false values
    return all(comb[i] == checks[i] for i in range(7))

# Generate and check all possible combinations
solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if check_combination(combination):
        solutions.append(combination)

print(len(solutions))