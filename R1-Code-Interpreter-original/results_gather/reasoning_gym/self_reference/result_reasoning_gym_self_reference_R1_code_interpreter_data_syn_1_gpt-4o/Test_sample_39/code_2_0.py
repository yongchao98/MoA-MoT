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
    # Check each condition
    cond1 = (T >= 4)
    cond2 = (F <= 2)
    cond3 = (T == 1)
    cond4 = (T == 0)
    cond5 = (T == 1 or T == 0) and not (T == 1 and T == 0)
    cond6 = is_prime(T)
    cond7 = is_composite(F)
    
    # Check if the combination of true and false statements satisfies the conditions
    if (cond1 and cond2 and cond3 and cond4 and cond5 and cond6 and cond7):
        possible_solutions += 1

print(possible_solutions)