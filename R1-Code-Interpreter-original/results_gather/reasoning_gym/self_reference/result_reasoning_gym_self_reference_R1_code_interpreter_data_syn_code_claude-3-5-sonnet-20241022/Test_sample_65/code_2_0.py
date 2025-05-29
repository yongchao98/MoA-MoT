def check_statements(true_count, statements):
    new_statements = statements.copy()
    
    # Statement 1: At least 5 statements are true
    new_statements[0] = (true_count >= 5)
    
    # Statement 2: At most 0 statements are false
    new_statements[1] = ((7 - true_count) <= 0)
    
    # Statement 3: Exactly 5 statements are true
    new_statements[2] = (true_count == 5)
    
    # Statement 4: Exactly 5 statements are false
    new_statements[3] = ((7 - true_count) == 5)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    new_statements[4] = (new_statements[2] != new_statements[3] and (new_statements[2] or new_statements[3]))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    new_statements[5] = (true_count in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}
    new_statements[6] = ((7 - true_count) in composites)
    
    return new_statements

def find_unique_solutions():
    solutions = set()
    for i in range(128):  # 2^7 possibilities
        statements = [False] * 7
        for j in range(7):
            statements[j] = bool((i >> j) & 1)
        
        initial_true_count = sum(statements)
        final_statements = check_statements(initial_true_count, statements)
        final_true_count = sum(final_statements)
        
        if final_true_count == initial_true_count:
            solution_tuple = tuple(final_statements)
            if solution_tuple not in solutions:
                solutions.add(solution_tuple)
                print(f"Unique solution found: {list(final_statements)} with {final_true_count} true statements")
    
    return len(solutions)

print(f"\nTotal number of unique solutions: {find_unique_solutions()}")