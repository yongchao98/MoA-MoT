def check_statements(true_statements):
    # Count true and false statements
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Initialize all statements as False
    is_valid = [False] * 7
    
    # Statement 1: At least 4 are true
    is_valid[0] = (num_true >= 4)
    
    # Statement 2: At most 1 is false
    is_valid[1] = (num_false <= 1)
    
    # Statement 3: Exactly 0 are true
    is_valid[2] = (num_true == 0)
    
    # Statement 4: Exactly 3 are false
    is_valid[3] = (num_false == 3)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    stmt3_true = 2 in true_statements
    stmt4_true = 3 in true_statements
    is_valid[4] = (stmt3_true != stmt4_true)  # XOR
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    is_valid[5] = (num_true in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    is_valid[6] = (num_false in composites)
    
    # Check if the combination is consistent
    return all(i in true_statements for i in range(7) if is_valid[i]) and \
           all(i not in true_statements for i in range(7) if not is_valid[i])

# Generate all possible combinations and check them
solutions = []
for i in range(2**7):
    true_statements = {j for j in range(7) if (i & (1 << j))}
    if check_statements(true_statements):
        solutions.append(true_statements)

print(len(solutions))