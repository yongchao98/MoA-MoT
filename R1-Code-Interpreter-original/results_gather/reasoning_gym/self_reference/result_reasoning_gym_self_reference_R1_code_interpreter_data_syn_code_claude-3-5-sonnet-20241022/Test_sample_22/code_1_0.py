def check_statements(true_count):
    false_count = 7 - true_count
    
    # Statement 1: At least 7 statements are true
    s1 = (true_count >= 7)
    
    # Statement 2: At most 2 statements are false
    s2 = (false_count <= 2)
    
    # Statement 3: Exactly 4 statements are true
    s3 = (true_count == 4)
    
    # Statement 4: Exactly 5 statements are false
    s4 = (false_count == 5)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    s5 = (s3 != s4) and (s3 or s4)
    
    # Statement 6: Number of true statements is prime
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    s6 = is_prime(true_count)
    
    # Statement 7: Number of false statements is composite
    def is_composite(n):
        return n > 1 and not is_prime(n)
    s7 = is_composite(false_count)
    
    statements = [s1, s2, s3, s4, s5, s6, s7]
    actual_true = sum(statements)
    
    if actual_true == true_count:
        print(f"Solution found with {true_count} true statements:")
        for i, s in enumerate(statements, 1):
            print(f"Statement {i}: {'True' if s else 'False'}")
        print()
        return True
    return False

solutions = sum(1 for i in range(8) if check_statements(i))
print(f"\nTotal number of solutions: {solutions}")