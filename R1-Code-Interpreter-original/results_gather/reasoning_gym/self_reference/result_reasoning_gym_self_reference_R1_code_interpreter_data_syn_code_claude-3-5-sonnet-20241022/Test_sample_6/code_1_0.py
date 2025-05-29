def verify_solution_detailed(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n**0.5) + 1):
            if n % i == 0: return False
        return True
    
    # Check each statement's truth value
    actual_values = [
        num_true >= 4,                    # Statement 1
        num_false <= 6,                   # Statement 2
        num_true == 4,                    # Statement 3
        num_false == 3,                   # Statement 4
        not (3 in true_statements and 4 in true_statements) and 
        (3 in true_statements or 4 in true_statements),  # Statement 5
        is_prime(num_true),               # Statement 6
        num_false in {4,6,8,9}           # Statement 7
    ]
    
    print(f"\nChecking solution {sorted(list(true_statements))}:")
    print(f"Number of true statements: {num_true}")
    print(f"Number of false statements: {num_false}")
    
    for i, value in enumerate(actual_values):
        print(f"Statement {i}: Should be {'true' if i in true_statements else 'false'}, Is {value}")
    
    return all((i in true_statements) == value for i, value in enumerate(actual_values))

solutions = [set(), {1,5}, {1,5,6}]
valid = sum(1 for s in solutions if verify_solution_detailed(s))
print(f"\nTotal confirmed valid solutions: {valid}")