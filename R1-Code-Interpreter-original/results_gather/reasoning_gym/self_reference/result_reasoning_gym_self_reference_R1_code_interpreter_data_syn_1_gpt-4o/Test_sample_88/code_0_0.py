def is_prime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

possible_solutions = 0

for true_statements in range(8):
    false_statements = 7 - true_statements
    
    # Check each statement
    s1 = true_statements >= 6
    s2 = false_statements <= 4
    s3 = true_statements == 5
    s4 = false_statements == 2
    s5 = (s3 or s4) and not (s3 and s4)
    s6 = is_prime(true_statements)
    s7 = is_composite(false_statements)
    
    # Count the number of true statements
    true_count = sum([s1, s2, s3, s4, s5, s6, s7])
    
    # Check if the number of true statements matches the condition
    if true_count == true_statements:
        possible_solutions += 1

print(possible_solutions)