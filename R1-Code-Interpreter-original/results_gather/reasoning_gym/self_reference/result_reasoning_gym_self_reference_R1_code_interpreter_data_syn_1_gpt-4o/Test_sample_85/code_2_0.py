def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements():
    possible_solutions = 0
    for num_true in range(8):
        num_false = 7 - num_true
        
        # Check each statement
        s1 = num_true >= 4
        s2 = num_false <= 5
        s3 = num_true == 4
        s4 = num_false == 7
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = is_prime(num_true)
        s7 = is_composite(num_false)
        
        # Determine which statements are true
        true_statements = [s1, s2, s3, s4, s5, s6, s7]
        
        # Count the number of true statements
        if true_statements.count(True) == num_true:
            possible_solutions += 1

    return possible_solutions

print(check_statements())