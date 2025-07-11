# The problem is to evaluate the infinite product from n=3 to infinity of (1 - z^3/n^3).
# The solution is derived using the product representation of the Gamma function.
# The final expression involves the Gamma function and the complex cube roots of unity.
# This script prints the final formula as a string.

# The numbers in the final equation are derived from the starting index 3,
# which leads to factoring out terms for n=1 and n=2.
# 1^3 = 1
# 2^3 = 8
# The cube roots of unity are represented using exp(2*pi*i/3) and exp(4*pi*i/3).
# The numbers in the final equation are 1, 8, 1, 2, 3, 4.

print("1 / ( (1 - z**3) * (1 - z**3/8) * Gamma(1 - z) * Gamma(1 - z*exp(2*pi*i/3)) * Gamma(1 - z*exp(4*pi*i/3)) )")