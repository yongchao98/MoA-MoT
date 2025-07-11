import math

def euler_totient(n):
    """Computes Euler's totient function."""
    if n <= 0:
        return 0
    # For a prime number p, phi(p) = p - 1.
    # Since 31 is a prime number, we can use this simple formula.
    # A more general implementation is provided for completeness, but not strictly necessary here.
    count = 0
    for i in range(1, n + 1):
        if math.gcd(n, i) == 1:
            count += 1
    return count

# The order of the cyclic group C
n = 31

# The problem reduces to calculating the order of Aut(C), which is phi(n).
# C is Z_31, so we calculate phi(31).
# Since 31 is a prime number, phi(31) = 31 - 1.
p = 31
one = 1
result = p - one

# The problem asks to output each number in the final equation.
# The final equation can be represented as 31 - 1 = 30.
print(f"The sum is calculated by the equation {p} - {one} = {result}.")

<<<30>>>