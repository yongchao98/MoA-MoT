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

# Check all possible numbers of true statements
for true_statements in range(8):
    false_statements = 7 - true_statements
    
    # Check if the configuration satisfies all conditions
    if (true_statements >= 1 and
        false_statements <= 3 and
        is_prime(true_statements) and
        is_composite(false_statements) and
        ((true_statements == 6 and false_statements == 1) or
         (true_statements == 4 and false_statements == 3)) and
        ((true_statements == 6) != (true_statements == 4))):
        possible_solutions += 1

print(possible_solutions)