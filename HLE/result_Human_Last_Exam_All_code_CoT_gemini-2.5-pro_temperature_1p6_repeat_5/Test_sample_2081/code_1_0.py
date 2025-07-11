import math
from decimal import Decimal, getcontext

# Set precision for Decimal calculations to handle large numbers accurately
getcontext().prec = 100

# Plan:
# 1. Define the value of T and use it to find e^T and e^(2T).
#    Given T = ln(10^34), e^T = 10^34.
# 2. Use the derived formula R^2 = (e^(2T) + e^T) / 2 to find R^2.
# 3. Calculate R by taking the square root of R^2.
# 4. Print the final equation with all numbers included, as requested.

# Step 1: Define e^T and e^(2T) using the Decimal type for high precision.
e_T = Decimal(10)**34
e_2T = Decimal(10)**68
two = Decimal(2)

# Step 2: Calculate R^2 using the derived formula.
R_squared = (e_2T + e_T) / two

# Step 3: Calculate R.
R = R_squared.sqrt()

# Step 4: Output the calculation, showing each number in the equation.
# We format the large numbers in scientific notation for readability.
e_2T_str = f"{e_2T:.1e}"
e_T_str = f"{e_T:.1e}"
R_str = f"{R:.8e}"

print("The set of initial values (x_0, y_0, z_0) for which a solution exists forms a sphere x_0^2 + y_0^2 + z_0^2 = R^2.")
print("The radius R is calculated from the formula: R = sqrt((e^(2T) + e^T) / 2)")
print("\nSubstituting the numerical values:")
print(f"e^T = {e_T_str}")
print(f"e^(2T) = {e_2T_str}")
print("\nThe final equation for R is:")
print(f"R = sqrt(({e_2T_str} + {e_T_str}) / {two}) = {R_str}")
