def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(statements):
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement's validity
    valid = [False] * 7
    
    # Statement 1: At least 7 statements are true
    valid[0] = (true_count >= 7) == statements[0]
    
    # Statement 2: At most 2 statements are false
    valid[1] = (false_count <= 2) == statements[1]
    
    # Statement 3: Exactly 0 statements are true
    valid[2] = (true_count == 0) == statements[2]
    
    # Statement 4: Exactly 4 statements are false
    valid[3] = (false_count == 4) == statements[3]
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid[4] = ((statements[2] != statements[3]) and (statements[2] or statements[3])) == statements[4]
    
    # Statement 6: Number of true statements is prime
    valid[5] = is_prime(true_count) == statements[5]
    
    # Statement 7: Number of false statements is composite
    valid[6] = is_composite(false_count) == statements[6]
    
    return all(valid)

# Try all possible combinations
solutions = []
for i in range(2**7):
    statements = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(statements):
        solutions.append(statements)

print(len(solutions))
if len(solutions) > 0:
    print("Example solution:", solutions[0])