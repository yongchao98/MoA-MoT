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

def check_statements():
    possible_solutions = 0
    for true_count in range(8):
        false_count = 7 - true_count
        statement_1 = true_count >= 1
        statement_2 = false_count <= 3
        statement_3 = true_count == 4
        statement_4 = false_count == 3
        statement_5 = (statement_3 or statement_4) and not (statement_3 and statement_4)
        statement_6 = is_prime(true_count)
        statement_7 = is_composite(false_count)
        
        if (statement_1 and statement_2 and statement_5 and statement_6 and statement_7):
            possible_solutions += 1
    
    print(possible_solutions)

check_statements()