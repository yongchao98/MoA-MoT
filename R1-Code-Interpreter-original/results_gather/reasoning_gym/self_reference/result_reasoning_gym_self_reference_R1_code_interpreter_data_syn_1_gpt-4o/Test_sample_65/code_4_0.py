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
    statement1 = (T >= 5)
    statement2 = (T == 7)
    statement3 = (T == 5)
    statement4 = (T == 2)
    statement5 = (statement3 != statement4)
    statement6 = is_prime(T)
    statement7 = is_composite(F)
    
    # We need to find a consistent set of true/false values
    if (statement1 or statement2 or statement3 or statement4 or statement5 or statement6 or statement7) and \
       (not statement1 or not statement2 or not statement3 or not statement4 or not statement5 or not statement6 or not statement7):
        possible_solutions += 1

print(possible_solutions)