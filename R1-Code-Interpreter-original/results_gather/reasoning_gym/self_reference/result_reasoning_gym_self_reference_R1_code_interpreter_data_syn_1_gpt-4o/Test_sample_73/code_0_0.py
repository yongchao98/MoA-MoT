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

for T in range(8):  # T can be from 0 to 7
    F = 7 - T
    if (T >= 5 and F <= 6 and
        ((T == 5) != (F == 2)) and
        is_prime(T) and
        is_composite(F)):
        possible_solutions += 1

print(possible_solutions)