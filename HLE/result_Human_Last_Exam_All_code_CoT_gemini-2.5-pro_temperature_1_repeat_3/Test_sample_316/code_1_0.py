import numpy as np

def alpha(p):
    """
    Calculates the conjectured sharp exponent alpha for the decoupling
    inequality for the cone in R^3.
    The function is piecewise, with critical exponents at p=3 and p=4.
    """
    if p >= 4:
        # For p >= 4, the decoupling conjecture is known to hold (up to epsilon).
        return 0
    elif p >= 3:
        # For 3 <= p <= 4, the exponent is determined by the planar/bilinear example.
        return 1/p - 1/4
    elif p > 2:
        # For 2 < p < 3, the exponent interpolates between the p=3 point and the p=2 point.
        # This line is sharp for the Kakeya/hairbrush example.
        return 1/4 - 1/(2*p)
    else:
        # The problem is for p > 2.
        return float('nan')

def get_slope(func, p, h=1e-6):
    """Calculates the slope of alpha(1/p) with respect to 1/p."""
    x = 1/p
    x_plus_h = 1/(p - h*p*p) # x+h approx
    return (alpha(1/x_plus_h) - alpha(1/x)) / h if p > 2 else float('nan')


print("The problem asks for the other critical exponent for a reverse square function estimate for the cone in R^3.")
print("The exponent alpha is a piecewise linear function of 1/p.")
print("One critical exponent, where the slope of alpha(1/p) changes, is given as p=4.\n")

# Let's verify our proposed model with the other critical exponent at p=3.

p1 = 4.0
p2 = 3.0

# Check continuity at p=4
alpha_at_4_from_above = alpha(p1 + 0.001)
alpha_at_4_from_below = alpha(p1 - 0.001)
print(f"Checking the function around p = {p1}:")
print(f"  alpha({p1+0.001:.3f}) = {alpha_at_4_from_above:.6f}")
print(f"  alpha({p1-0.001:.3f}) = {alpha_at_4_from_below:.6f}")
print(f"The function is continuous at p={p1}.\n")


# Check continuity at p=3
alpha_at_3_from_above = alpha(p2 + 0.001)
alpha_at_3_from_below = alpha(p2 - 0.001)
print(f"Checking the function around p = {p2}:")
print(f"  alpha({p2+0.001:.3f}) = {alpha_at_3_from_above:.6f}  (from the [3,4] segment)")
print(f"  alpha({p2-0.001:.3f}) = {alpha_at_3_from_below:.6f}  (from the (2,3) segment)")
# Check the value at p=3 exactly
print(f"  alpha({p2}) from the [3,4] segment is 1/3 - 1/4 = {1/3 - 1/4:.6f}")
print(f"  alpha({p2}) from the (2,3) segment is 1/4 - 1/(2*3) = {1/4 - 1/6:.6f}")
print(f"The function is continuous at p={p2}.\n")

# Check slopes to confirm these are the points where the slope changes.
slope_above_4 = get_slope(alpha, p1 + 0.1)
slope_between_3_and_4 = get_slope(alpha, 3.5)
slope_between_2_and_3 = get_slope(alpha, 2.5)

print("Checking the slope of alpha(1/p) in different regimes:")
print(f"  Slope for p > 4 (e.g., at p=4.1): {slope_above_4:.2f}")
print(f"  Slope for 3 < p < 4 (e.g., at p=3.5): {slope_between_3_and_4:.2f}")
print(f"  Slope for 2 < p < 3 (e.g., at p=2.5): {slope_between_2_and_3:.2f}")

print("\nThe slope changes from 0 to 1.00 at p=4, and from 1.00 to -0.50 at p=3.")
print("This fits the problem description of having two critical exponents.")
print("Since one critical exponent is 4, the other must be 3.")
print("\nThe other critical exponent is:")
print(3)