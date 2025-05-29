def check_statements(truths):
    # Count number of true and false statements
    num_true = sum(truths)
    num_false = 7 - num_true
    
    # Check each statement's consistency
    # Statement 1: At least 7 statements are true
    if truths[0] and num_true < 7:
        return False
    if not truths[0] and num_true >= 7:
        return False
    
    # Statement 2: At most 5 statements are false
    if truths[1] and num_false > 5:
        return False
    if not truths[1] and num_false <= 5:
        return False
    
    # Statement 3: Exactly 3 statements are true
    if truths[2] and num_true != 3:
        return False
    if not truths[2] and num_true == 3:
        return False
    
    # Statement 4: Exactly 7 statements are false
    if truths[3] and num_false != 7:
        return False
    if not truths[3] and num_false == 7:
        return False
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    if truths[4] and not (truths[2] ^ truths[3]):
        return False
    if not truths[4] and (truths[2] ^ truths[3]):
        return False
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    if truths[5] and num_true not in primes:
        return False
    if not truths[5] and num_true in primes:
        return False
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    if truths[6] and num_false not in composites:
        return False
    if not truths[6] and num_false in composites:
        return False
    
    return True

# Generate all possible combinations and show valid solutions
solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    truths = [bool(i & (1 << j)) for j in range(7)]
    if check_statements(truths):
        print(f"Solution found: {truths}")
        print(f"Number of true statements: {sum(truths)}")
        print(f"Number of false statements: {7-sum(truths)}")
        print("---")
        solutions.append(truths)

print(f"Total number of solutions: {len(solutions)}")