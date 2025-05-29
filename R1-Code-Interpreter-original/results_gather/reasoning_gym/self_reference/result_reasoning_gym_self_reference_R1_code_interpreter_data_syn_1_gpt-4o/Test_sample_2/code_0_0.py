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
    if (T >= 4 and F <= 1 and (T == 5 or F == 2) and (T == 5) and is_prime(T) and is_composite(F)):
        possible_solutions += 1

print(possible_solutions)