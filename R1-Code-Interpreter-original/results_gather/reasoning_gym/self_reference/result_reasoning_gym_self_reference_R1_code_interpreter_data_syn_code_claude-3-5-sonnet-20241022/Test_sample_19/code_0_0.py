def check_statements(statements):
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Statement 1: At least 7 statements are true
    if statements[0] and true_count < 7:
        return False
    if not statements[0] and true_count >= 7:
        return False
    
    # Statement 2: At most 1 statement is false
    if statements[1] and false_count > 1:
        return False
    if not statements[1] and false_count <= 1:
        return False
    
    # Statement 3: Exactly 1 statement is true
    if statements[2] and true_count != 1:
        return False
    if not statements[2] and true_count == 1:
        return False
    
    # Statement 4: Exactly 7 statements are false
    if statements[3] and false_count != 7:
        return False
    if not statements[3] and false_count == 7:
        return False
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    if statements[4] and not (statements[2] ^ statements[3]):
        return False
    if not statements[4] and (statements[2] ^ statements[3]):
        return False
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    if statements[5] and true_count not in primes:
        return False
    if not statements[5] and true_count in primes:
        return False
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    if statements[6] and false_count not in composites:
        return False
    if not statements[6] and false_count in composites:
        return False
    
    return True

# Try all possible combinations
solutions = []
for i in range(2**7):
    statements = [bool((i >> j) & 1) for j in range(7)]
    if check_statements(statements):
        solutions.append(statements)

print(f"Number of solutions: {len(solutions)}")