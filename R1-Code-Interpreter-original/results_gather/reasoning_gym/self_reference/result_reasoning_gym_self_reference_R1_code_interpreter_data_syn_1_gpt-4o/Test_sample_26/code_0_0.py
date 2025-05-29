from itertools import combinations

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

def count_valid_solutions():
    valid_solutions = 0
    for true_statements in range(1, 8):
        false_statements = 7 - true_statements
        if true_statements >= 4 and is_prime(true_statements) and is_composite(false_statements):
            if (true_statements == 5 and false_statements == 2) or (true_statements == 3 and false_statements == 4):
                valid_solutions += 1
    return valid_solutions

print(count_valid_solutions())