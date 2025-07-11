import math

# The R value where the standard logistic map is chaotic
R = 3.57

# The modified logistic map is X_n+1 = R * X_n * (1 - X_n) + (R-1)/2
# To find the equilibrium point X_eq, we solve the equation:
# X_eq = R * X_eq * (1 - X_eq) + (R-1)/2
# This rearranges to a quadratic equation:
# R*(X_eq)^2 - (R-1)*X_eq - (R-1)/2 = 0

# Coefficients for the quadratic formula a*x^2 + b*x + c = 0
a = R
b = -(R - 1)
c = -(R - 1) / 2

# Calculate the discriminant
discriminant = b**2 - 4 * a * c

# The equilibrium point must be real, so the discriminant must be non-negative
if discriminant >= 0:
    # We choose the positive root as it corresponds to the desired value
    equilibrium_point = (-b + math.sqrt(discriminant)) / (2 * a)

    print("The modified logistic map is given by the equation:")
    # Using f-string to output each number in the equation
    term1_val = R
    term2_val = R - 1
    term3_val = 2
    
    print(f"X_n+1 = {term1_val} * X_n * (1 - X_n) + ({term2_val}) / {term3_val}")
    print("\nFor R = 3.57, this equation replaces chaos with a stable equilibrium point.")
    print(f"The calculated equilibrium point is approximately: {equilibrium_point}")
else:
    print("No real equilibrium point exists for the given R.")
