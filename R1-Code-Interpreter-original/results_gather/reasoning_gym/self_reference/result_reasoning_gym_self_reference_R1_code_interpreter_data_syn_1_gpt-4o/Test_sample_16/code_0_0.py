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
    
    # Check conditions
    if true_statements >= 6:  # Statement 1
        if true_statements == 5:  # Statement 3
            if false_statements == 2:  # Statement 4
                if is_prime(true_statements):  # Statement 6
                    if is_composite(false_statements):  # Statement 7
                        possible_solutions += 1

print(possible_solutions)