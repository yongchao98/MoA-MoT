def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

possible_solutions = 0

for T in range(8):
    F = 7 - T
    if (T >= 3 and F <= 3 and
        ((T == 4 and F != 3) or (F == 3 and T != 4)) and
        is_prime(T) and is_composite(F)):
        possible_solutions += 1

print(possible_solutions)