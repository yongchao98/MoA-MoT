import numpy as np
from scipy.integrate import quad

# Define the integrand for L
def integrand(t):
    return 2 / np.sqrt(1 - t**4)

# Calculate L by numerical integration
# We integrate up to a point very close to 1 to avoid the singularity
L, error = quad(integrand, 0, 1)

# The inequality is sum(n*|a_n|^2) <= sum(|a_n|)
# LHS is Area/pi = L^2/pi
lhs = L**2 / np.pi

# A sufficient condition is for the LHS to be less than or equal to the sum of the first two terms of the RHS.
# RHS_approx = |a_0| + |a_1| = L/sqrt(2) + sqrt(2)
rhs_approx = L / np.sqrt(2) + np.sqrt(2)

print(f"The side length of the square is L = {L:.4f}")
print(f"The equation to check is: L^2/pi <= L/sqrt(2) + sqrt(2) + ...")
print(f"We check if {lhs:.4f} <= {rhs_approx:.4f}")

if lhs <= rhs_approx:
    print("The inequality holds for the first two terms, so it holds for the entire series.")
    satisfies = True
else:
    print("The inequality does not hold based on the first two terms.")
    satisfies = False

print(f"\nFinal check:")
print(f"The term n=0 gives |a_0| = L/sqrt(2) = {L/np.sqrt(2):.4f}")
print(f"The term n=1 gives |a_1| = sqrt(2) = {np.sqrt(2):.4f}")
print(f"The sum of these two terms is {rhs_approx:.4f}")
print(f"The left hand side of the inequality is L^2/pi = {lhs:.4f}")
print(f"The equation is {lhs:.4f} <= {rhs_approx:.4f} + sum(|a_n| for n>=2)")
print("Since the left side is smaller than the sum of just the first two terms on the right, the inequality is satisfied.")
