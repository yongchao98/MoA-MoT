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

for true_count in range(1, 8):
    false_count = 7 - true_count
    
    # Check each statement
    statement_1 = true_count >= 4
    statement_2 = false_count <= 6
    statement_3 = true_count == 4
    statement_4 = false_count == 3
    statement_5 = (statement_3 or statement_4) and not (statement_3 and statement_4)
    statement_6 = is_prime(true_count)
    statement_7 = is_composite(false_count)
    
    # Count the number of true statements
    true_statements = sum([statement_1, statement_2, statement_3, statement_4, statement_5, statement_6, statement_7])
    
    # Check if the number of true statements matches the true_count
    if true_statements == true_count:
        possible_solutions += 1

print(possible_solutions)