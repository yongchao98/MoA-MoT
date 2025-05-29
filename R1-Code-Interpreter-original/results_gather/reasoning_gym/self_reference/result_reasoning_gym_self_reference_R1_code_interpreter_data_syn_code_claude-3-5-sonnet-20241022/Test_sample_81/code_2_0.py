def verify_solution(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # For readability, let's check each statement
    s1 = "At least 7 statements are true"
    s2 = "At most 5 statements are false"
    s3 = "Exactly 3 statements are true"
    s4 = "Exactly 0 statements are false"
    s5 = "Either Statement 3 or Statement 4 is true, but not both"
    s6 = "The number of true statements is a prime number"
    s7 = "The number of false statements is a composite number"
    
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n**0.5) + 1):
            if n % i == 0: return False
        return True
    
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    # Check each solution
    statements = [i in true_statements for i in range(7)]
    checks = [
        num_true >= 7,
        num_false <= 5,
        num_true == 3,
        num_false == 0,
        (3 in true_statements) != (4 in true_statements),
        is_prime(num_true),
        is_composite(num_false)
    ]
    
    print(f"\nAnalyzing solution with true statements: {[x+1 for x in sorted(list(true_statements))]}")
    print(f"Number of true statements: {num_true}")
    print(f"Number of false statements: {7-num_true}")
    for i, (stmt, check, is_true) in enumerate(zip(
        [s1,s2,s3,s4,s5,s6,s7], 
        checks, 
        statements
    )):
        print(f"Statement {i+1}: {'True' if is_true else 'False'}, Claim is {'True' if check else 'False'}")
        if is_true != check:
            print(f"INCONSISTENCY found in statement {i+1}")

# Test each solution
solutions = [set(), {1,5}, {6}]
for sol in solutions:
    verify_solution(sol)