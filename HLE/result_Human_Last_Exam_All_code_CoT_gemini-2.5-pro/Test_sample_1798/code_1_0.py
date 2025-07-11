import math

# Step 1: Define the circle's parameters based on the derived equation
# Circle equation: (x - h)^2 + (y - k)^2 = r2
h = 5.5
k = -1.5
r2 = 22.5
r = math.sqrt(r2)

# The function is f(x) = k + sqrt(r2 - (x-h)^2)
# The PDF is d_X(x) = alpha * f(x)

# Step 2: Define a helper function to compute the integral of the circular part
# This function calculates the indefinite integral of sqrt(r^2 - u^2)
def G(u, r_val, r2_val):
    """Computes the indefinite integral of sqrt(r^2 - u^2)"""
    # Clamp the argument of asin to avoid math domain errors due to precision
    asin_arg = max(-1.0, min(1.0, u / r_val))
    return (u / 2.0) * math.sqrt(r2_val - u**2) + (r2_val / 2.0) * math.asin(asin_arg)

# Step 3: Calculate the total integral of f(x) from 1 to 10 to find alpha
# I_total = integral from 1 to 10 of (k + sqrt(r2 - (x-h)^2)) dx
# I_total = k*(10-1) + integral from 1 to 10 of sqrt(r2 - (x-h)^2) dx
# Let u = x - h. When x=1, u=1-5.5=-4.5. When x=10, u=10-5.5=4.5
integral_sqrt_total = G(10.0 - h, r, r2) - G(1.0 - h, r, r2)
I_total = k * (10.0 - 1.0) + integral_sqrt_total

# Calculate alpha
alpha = 1.0 / I_total

# Step 4: Calculate the partial integral of f(x) from 1 to 3 for the probability
# P(X<3) = alpha * integral from 1 to 3 of f(x) dx
# I_partial = integral from 1 to 3 of (k + sqrt(r2 - (x-h)^2)) dx
# I_partial = k*(3-1) + integral from 1 to 3 of sqrt(r2 - (x-h)^2) dx
# Let u = x - h. When x=1, u=-4.5. When x=3, u=3-5.5=-2.5
integral_sqrt_partial = G(3.0 - h, r, r2) - G(1.0 - h, r, r2)
I_partial = k * (3.0 - 1.0) + integral_sqrt_partial

# Calculate the final probability
prob_X_less_than_3 = alpha * I_partial

# Step 5: Output the results
print(f"The equation of the circle is (x - {h})^2 + (y - ({k}))^2 = {r2}")
print(f"The probability density function is d_X(x) = alpha * f(x), where f(x) = {k} + sqrt({r2} - (x - {h})^2)")
print("\n--- Calculated Values ---")
print(f"The value of the normalization constant alpha is: {alpha}")
print(f"The probability P(X < 3) is calculated as alpha * (integral of f(x) from 1 to 3).")
print(f"P(X < 3) = {alpha} * {I_partial}")
print(f"The final value for P(X < 3) is: {prob_X_less_than_3}")
