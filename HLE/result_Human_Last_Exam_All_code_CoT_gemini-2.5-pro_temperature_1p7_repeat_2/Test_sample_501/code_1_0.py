# This script presents the derived force law for a thermally isolated polymer chain.

# The problem asks for the force of attraction, F, between the polymer ends.
# This force is a function of several parameters:
# x: the separation of the ends
# n: the number of segments in the polymer chain
# l: the length of a single segment
# E(0): the kinetic energy of the polymer at zero extension

# Based on statistical mechanics in the microcanonical ensemble, the derived force law is:
# F(x) = - (2 * E(0) * x) / (n^2 * l^2)
# The negative sign indicates that it is an attractive (restoring) force.

# The code below will print this formula, highlighting the numerical coefficients.
coefficient_numerator = 2
power_of_n = 2
power_of_l = 2

print("The force law for the attraction between the polymer ends is F(x), where:")
print(f"F(x) = - ({coefficient_numerator} * E(0) * x) / (n^{power_of_n} * l^{power_of_l})")
print("\nThis law is valid for small extensions x and large n, and defines the following variables:")
print(" F(x): The force of attraction between the polymer ends.")
print(" E(0): The kinetic energy of the polymer when the extension x is zero.")
print(" x: The separation of the polymer ends.")
print(" n: The number of segments in the chain.")
print(" l: The length of each segment (using 'l' to represent Greek letter ell).")
