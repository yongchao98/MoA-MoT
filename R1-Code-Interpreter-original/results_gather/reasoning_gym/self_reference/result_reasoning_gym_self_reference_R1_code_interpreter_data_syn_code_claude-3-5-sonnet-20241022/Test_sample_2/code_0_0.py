def check_statements(true_count):
    # Returns if statements are true based on number of true statements
    statements = [False] * 7
    
    # Statement 1: At least 4 are true
    statements[0] = (true_count >= 4)
    
    # Statement 2: At most 1 is false
    statements[1] = ((7 - true_count) <= 1)
    
    # Statement 3: Exactly 5 are true
    statements[2] = (true_count == 5)
    
    # Statement 4: Exactly 2 are false
    statements[3] = ((7 - true_count) == 2)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    statements[4] = (statements[2] != statements[3]) and (statements[2] or statements[3])
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    statements[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for false statements (0-7)
    statements[6] = ((7 - true_count) in composites)
    
    # Count true statements in our evaluation
    actual_true_count = sum(statements)
    
    # Solution is valid if the number of true statements matches our assumption
    return actual_true_count == true_count

# Check all possible numbers of true statements (0 to 7)
solutions = []
for true_count in range(8):
    if check_statements(true_count):
        solutions.append(true_count)

print(f"Number of solutions: {len(solutions)}")
print(f"Solutions (number of true statements): {solutions}")