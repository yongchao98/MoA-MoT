def check_statements(true_count, statements):
    # Statement 1: At least 5 statements are true
    statements[0] = (true_count >= 5)
    
    # Statement 2: At most 0 statements are false
    statements[1] = ((7 - true_count) <= 0)
    
    # Statement 3: Exactly 5 statements are true
    statements[2] = (true_count == 5)
    
    # Statement 4: Exactly 5 statements are false
    statements[3] = ((7 - true_count) == 5)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    statements[4] = (statements[2] != statements[3] and (statements[2] or statements[3]))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    statements[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    statements[6] = ((7 - true_count) in composites)
    
    # Count true statements after evaluation
    return sum(1 for s in statements if s)

def find_solutions():
    solutions = []
    # Try all possible initial combinations
    for i in range(128):  # 2^7 possibilities
        statements = [False] * 7
        # Convert number to binary to get initial true/false values
        for j in range(7):
            statements[j] = bool((i >> j) & 1)
        
        initial_true_count = sum(1 for s in statements if s)
        final_true_count = check_statements(initial_true_count, statements)
        
        # If the number of true statements matches before and after evaluation
        if final_true_count == initial_true_count:
            solutions.append(statements[:])
    
    return len(solutions)

print(find_solutions())