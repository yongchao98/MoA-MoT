# Check all possible values of T (number of true statements) and F (number of false statements)
# to find consistent solutions.

def is_prime(n):
    """Check if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    """Check if a number is composite."""
    return n > 1 and not is_prime(n)

possible_solutions = 0

for T in range(8):  # T can be from 0 to 7
    F = 7 - T
    # Check the conditions
    if (T >= 3 and F <= 5 and
        ((T == 7) != (T == 6)) and  # Statement 5: Either Statement 3 or 4 is true, but not both
        is_prime(T) and
        is_composite(F)):
        possible_solutions += 1

print(possible_solutions)