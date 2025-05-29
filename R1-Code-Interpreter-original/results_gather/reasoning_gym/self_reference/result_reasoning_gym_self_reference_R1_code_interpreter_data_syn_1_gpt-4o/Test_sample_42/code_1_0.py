def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

# Check the specific case where T = 5 and F = 2
T = 5
F = 7 - T

# Verify all conditions
if (T >= 3 and F <= 3 and
    (T == 7 or T == 5) and
    is_prime(T) and
    not is_composite(F)):
    solutions = 1
else:
    solutions = 0

print(solutions)