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

for true_count in range(8):
    false_count = 7 - true_count
    
    # Evaluate each statement
    s1 = (true_count >= 7)
    s2 = (false_count <= 2)
    s3 = (true_count == 4)
    s4 = (false_count == 5)
    s5 = (s3 != s4)
    s6 = is_prime(true_count)
    s7 = is_composite(false_count)
    
    # Count the number of true statements
    true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
    
    # Check if the number of true statements matches the true_count
    if true_statements == true_count:
        possible_solutions += 1

print(possible_solutions)