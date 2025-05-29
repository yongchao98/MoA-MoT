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
    if (T >= 4 and T <= 7) and (F <= 2) and (T == 1) and (T == 0) and (T == 1 or T == 0) and (is_prime(T)) and (is_composite(F)):
        possible_solutions += 1

print(possible_solutions)