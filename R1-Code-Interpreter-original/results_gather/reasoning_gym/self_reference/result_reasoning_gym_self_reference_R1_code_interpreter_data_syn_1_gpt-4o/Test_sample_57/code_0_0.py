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
    
    # Check each condition
    cond1 = true_statements >= 6
    cond2 = false_statements <= 5
    cond3 = true_statements == 3
    cond4 = true_statements == 6
    cond5 = (cond3 or cond4) and not (cond3 and cond4)
    cond6 = is_prime(true_statements)
    cond7 = is_composite(false_statements)
    
    # All conditions must be satisfied
    if cond1 and cond2 and cond5 and cond6 and cond7:
        possible_solutions += 1

print(possible_solutions)