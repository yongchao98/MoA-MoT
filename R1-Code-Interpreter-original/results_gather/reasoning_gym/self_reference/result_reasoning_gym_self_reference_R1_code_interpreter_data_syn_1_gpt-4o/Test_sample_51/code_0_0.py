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

# Check scenario where Statement 3 is true
true_statements = 6
false_statements = 1

if (true_statements >= 1 and
    false_statements <= 3 and
    is_prime(true_statements) and
    is_composite(false_statements)):
    possible_solutions += 1

# Check scenario where Statement 4 is true
true_statements = 4
false_statements = 3

if (true_statements >= 1 and
    false_statements <= 3 and
    is_prime(true_statements) and
    is_composite(false_statements)):
    possible_solutions += 1

print(possible_solutions)