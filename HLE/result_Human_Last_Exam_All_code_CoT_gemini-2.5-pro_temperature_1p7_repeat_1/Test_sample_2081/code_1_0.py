import math

# Step 1: Define the constants based on the problem statement.
# T = ln(10^34), so e^T = 10^34.
eT = 10**34
e2T = eT**2

# Step 2: Calculate the square of the radius, R^2.
# R^2 = 0.5 * (e^(2T) + e^T)
# In Python, integer arithmetic is exact, so the addition is performed with full precision
# before being converted to a float for the multiplication.
# However, due to the limits of standard float precision, the term eT is negligible
# compared to e2T. float(e2T + eT) will be the same as float(e2T).
R_squared = 0.5 * (e2T + eT)

# Step 3: Calculate the radius R by taking the square root.
R = math.sqrt(R_squared)

# Step 4: Print the final equation with the computed values.
# The instruction is to "output each number in the final equation".
print(f"The equation for the radius R is: R = sqrt(0.5 * (e^(2T) + e^T))")
print(f"where:")
# Use scientific notation for clarity.
print(f"e^T = {eT:.0e}")
print(f"e^(2T) = {e2T:.0e}")

# The addition (e2T + eT) is handled by Python's arbitrary-precision integers first.
# The result is then converted to a float.
final_val_in_sqrt = e2T + eT
print(f"R = sqrt(0.5 * ({final_val_in_sqrt:.0e}))")
print(f"R = {R}")