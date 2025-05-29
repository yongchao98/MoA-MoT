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
    statement1 = T >= 3
    statement2 = F <= 7
    statement3 = T == 2
    statement4 = F == 2
    statement5 = (statement3 or statement4) and not (statement3 and statement4)
    statement6 = is_prime(T)
    statement7 = is_composite(F)
    
    if statement1 and statement2 and statement5 and statement6 and statement7:
        possible_solutions += 1

print(possible_solutions)