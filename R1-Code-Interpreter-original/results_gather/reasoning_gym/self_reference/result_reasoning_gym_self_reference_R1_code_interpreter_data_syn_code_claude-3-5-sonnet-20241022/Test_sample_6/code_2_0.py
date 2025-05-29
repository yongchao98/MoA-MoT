def verify_and_explain(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n**0.5) + 1):
            if n % i == 0: return False
        return True
    
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    # Evaluate each statement
    evaluations = [
        (num_true >= 4, "At least 4 are true"),
        (num_false <= 6, "At most 6 are false"),
        (num_true == 4, "Exactly 4 are true"),
        (num_false == 3, "Exactly 3 are false"),
        (((2 in true_statements) != (3 in true_statements)), "Either stmt3 or stmt4 is true but not both"),
        (is_prime(num_true), "Number of true statements is prime"),
        (is_composite(num_false), "Number of false statements is composite")
    ]
    
    print(f"\nAnalyzing solution with {num_true} true statements: {sorted(list(true_statements))}")
    for i, (eval_result, description) in enumerate(evaluations):
        status = "true" if i in true_statements else "false"
        actual = "true" if eval_result else "false"
        print(f"Statement {i}: Should be {status}, Actually is {actual} - {description}")
    
    return all((i in true_statements) == eval_result for i, (eval_result, _) in enumerate(evaluations))

solutions = [set(), {0,1,2,3}, {1,5}, {1,5,6}]
valid = sum(1 for s in solutions if verify_and_explain(s))
print(f"\nTotal valid solutions: {valid}")