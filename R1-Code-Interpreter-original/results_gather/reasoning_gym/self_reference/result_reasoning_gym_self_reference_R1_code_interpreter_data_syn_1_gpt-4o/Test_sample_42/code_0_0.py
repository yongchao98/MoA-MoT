def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

# Possible number of true statements
possible_T = [3, 4, 5, 6, 7]

# Check each possible number of true statements
solutions = 0
for T in possible_T:
    F = 7 - T
    if (T >= 3 and F <= 3 and
        (T == 7 or T == 5) and
        is_prime(T) and
        is_composite(F)):
        solutions += 1

print(solutions)