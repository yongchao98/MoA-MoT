import numpy as np

def analyze_torus_mean_curvature(R, a):
    """
    Analyzes the mean curvature of a torus of revolution.
    
    A torus can be parameterized by major radius R and minor radius a.
    Its mean curvature H (with respect to the inward normal) is given by the formula:
    H(u) = (R + 2*a*cos(u)) / (2*a*(R + a*cos(u)))
    where u is the angle for the circular cross-section.

    The condition that the mean curvature never vanishes is H(u) != 0 for all u.
    """
    
    print(f"Analyzing a torus with major radius R = {R} and minor radius a = {a}")

    # Condition for the torus to be well-defined (no self-intersection)
    if not R > a > 0:
        print("Invalid parameters: we must have R > a > 0.")
        return

    # The denominator of H(u) is denom = 2*a*(R + a*cos(u)).
    # Since R > a, the term R + a*cos(u) is always positive.
    # Its minimum value is R - a.
    denom_min = 2 * a * (R - a)
    print(f"The denominator is always positive, with a minimum value of {denom_min:.4f}")

    # The sign of H(u) is determined by the numerator, N(u) = R + 2*a*cos(u).
    # H(u) is zero only if N(u) is zero.
    # We need to check if N(u) can be zero.
    # The minimum value of N(u) occurs at u = pi (cos(u) = -1).
    numerator_min = R + 2 * a * (-1)
    
    # The maximum value of N(u) occurs at u = 0 (cos(u) = 1).
    numerator_max = R + 2 * a * (1)

    print("\nThe mean curvature H(u) is proportional to the numerator N(u) = R + 2*a*cos(u).")
    print("The final equation for the numerator is N(u) = {} + {}*cos(u)".format(R, 2*a))
    print(f"The minimum value of the numerator is N_min = {R} - {2*a} = {numerator_min:.4f}")
    print(f"The maximum value of the numerator is N_max = {R} + {2*a} = {numerator_max:.4f}")

    # If we choose R > 2*a, the minimum value of the numerator will be positive.
    if R > 2 * a:
        print("\nSince R > 2*a, the minimum value of the numerator is positive.")
        print("Therefore, the numerator is always positive, and the mean curvature H(u) never vanishes.")
        print("This shows that a surface of genus 1 can have non-vanishing mean curvature.")
    else:
        print("\nSince R <= 2*a, the numerator can be zero or change sign.")
        print("To ensure non-vanishing mean curvature, one must choose R > 2*a.")

# Let's choose specific parameters that satisfy the condition R > 2a.
# This choice demonstrates that a genus-1 surface can have mean curvature that never vanishes.
R_val = 3.0
a_val = 1.0
analyze_torus_mean_curvature(R_val, a_val)
