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

for T in range(8):
    F = 7 - T
    # Check the conditions for each statement
    statement_1 = (T >= 5)
    statement_2 = (F == 0)
    statement_3 = (T == 5)
    statement_4 = (F == 4)
    statement_5 = (statement_3 != statement_4)
    statement_6 = is_prime(T)
    statement_7 = is_composite(F)
    
    # Check if the combination of true/false statements is consistent
    if (statement_1 + statement_2 + statement_3 + statement_4 + statement_5 + statement_6 + statement_7) == T:
        possible_solutions += 1

print(possible_solutions)