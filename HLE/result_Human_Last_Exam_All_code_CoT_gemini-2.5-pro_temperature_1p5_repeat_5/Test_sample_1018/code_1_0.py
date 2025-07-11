# Define the genera of the two surfaces.
g1 = 31
g2 = 17

# The problem is to compute the simplicial volume of the product space Sigma_g1 x Sigma_g2.
# Let M = Sigma_g1 and N = Sigma_g2.

# A fundamental theorem by Gromov states that if a closed aspherical manifold X
# has a fundamental group pi_1(X) containing a subgroup isomorphic to Z^2,
# then its simplicial volume ||X|| is 0.

# 1. Asphericity: A surface Sigma_g is aspherical for g >= 1.
#    Since g1=31 and g2=17 are both >= 1, the product manifold M x N is aspherical.

# 2. Fundamental Group: The fundamental group of the product is the product of the fundamental groups:
#    pi_1(M x N) = pi_1(M) x pi_1(N).
#    Since g1 >= 1 and g2 >= 1, we can pick non-trivial elements a in pi_1(M) and b in pi_1(N).
#    The elements (a, id) and (id, b) commute and generate a subgroup isomorphic to Z^2 in pi_1(M x N).

# Based on this theorem, the simplicial volume of the product space is 0.
result = 0

# Print the final equation as requested, showing the numbers involved.
print(f"The simplicial volume of the product of surfaces \u03A3_{g1} and \u03A3_{g2} is given by the equation:")
print(f"||\u03A3_{g1} \u00D7 \u03A3_{g2}|| = {result}")

# Per the instruction "output each number in the final equation!", we print the numbers g1, g2, and the result.
print("\nNumbers in the final equation:")
print(f"Genus 1: {g1}")
print(f"Genus 2: {g2}")
print(f"Result: {result}")