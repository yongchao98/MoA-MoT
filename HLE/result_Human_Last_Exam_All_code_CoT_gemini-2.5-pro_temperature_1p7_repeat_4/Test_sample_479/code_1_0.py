import numpy as np

def check_torus_mean_curvature(R, r):
    """
    Checks if a torus with major radius R and minor radius r has positive mean curvature everywhere.

    The mean curvature H of a torus of revolution has the same sign as the expression
    N(u) = R + 2*r*cos(u), where u is the poloidal angle.
    The denominator in the formula for H is always positive for a valid torus (R > r).
    For H to be positive everywhere, N(u) must be positive for all u.
    The minimum value of N(u) occurs when cos(u) = -1, which gives N_min = R - 2*r.
    So, H > 0 everywhere if and only if R - 2*r > 0, or R > 2*r.
    """
    print(f"Checking a torus with major radius R = {R} and minor radius r = {r}")
    print("The mean curvature H has the same sign as N(u) = R + 2*r*cos(u).")
    
    # Calculate the minimum value of the numerator N(u)
    min_numerator_val = R - 2 * r
    
    print(f"The final equation to check is if the minimum value of the numerator is positive.")
    print(f"Minimum numerator value = R - 2*r = {R} - 2*{r} = {min_numerator_val}")

    if min_numerator_val > 0:
        print(f"Since the minimum value {min_numerator_val} is positive, the mean curvature H is always positive.")
        print("This demonstrates that a surface of genus 1 (a torus) can satisfy the given condition.")
    else:
        print(f"Since the minimum value {min_numerator_val} is not positive, the mean curvature H is not always positive for this torus.")
        print("This example does not satisfy the condition.")

# We need to find a torus where R > 2r to satisfy the condition. Let's pick R=3, r=1.
check_torus_mean_curvature(3, 1)

print("\nBased on the demonstration that genus 1 is possible, and a theorem from differential geometry stating that the genus cannot exceed 1 for a surface with positive mean curvature in R^3, the maximal genus is 1.")
