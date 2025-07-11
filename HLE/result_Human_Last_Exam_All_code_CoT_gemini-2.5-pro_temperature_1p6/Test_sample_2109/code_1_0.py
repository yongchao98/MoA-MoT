import sympy as sp

# Step 1 & 2: Define the heat flux function and the expression for total energy
# From the derivation, Phi(z) = z/2 + 1.
# E_total = integral(Phi(Z(x)), (x, 0, 1)) = 1/2 * integral(Z(x), (x, 0, 1)) + 1
# where Z(x) = T1(sqrt(2)*x) + T2(x)

# Step 3: Use derived boundary values for T1 and T2
T1_0 = 0
T1_sqrt2 = sp.Rational(1, 3)

# Case A for T2
T2_0_A = 4
T2_1_A = 2

# Case B for T2
T2_0_B = sp.Rational(1, 3)
T2_1_B = sp.Rational(-5, 3)

# Step 4: Calculate Z(0) and Z(1) for both cases
# Z(0) = T1(0) + T2(0)
# Z(1) = T1(sqrt(2)) + T2(1)

Z0_A = T1_0 + T2_0_A
Z1_A = T1_sqrt2 + T2_1_A

Z0_B = T1_0 + T2_0_B
Z1_B = T1_sqrt2 + T2_1_B

# Approximate the integral of Z(x) using the trapezoidal rule
# integral(Z(x)) from 0 to 1 is approx (Z(0) + Z(1))/2
int_Z_A = (Z0_A + Z1_A) / 2
int_Z_B = (Z0_B + Z1_B) / 2

# Calculate the total energy for both cases
E_A = sp.Rational(1, 2) * int_Z_A + 1
E_B = sp.Rational(1, 2) * int_Z_B + 1

# The minimum energy is the smaller of the two values
min_E = min(E_A, E_B)

# Output the final calculation step by step
print("--- Calculation of Total Heat Energy ---")
print("\nCase A:")
print(f"Z(0) = T1(0) + T2(0) = {T1_0} + {T2_0_A} = {Z0_A}")
print(f"Z(1) = T1(sqrt(2)) + T2(1) = {T1_sqrt2} + {T2_1_A} = {Z1_A}")
print(f"Approximated Integral of Z(x) = (Z(0) + Z(1)) / 2 = ({Z0_A} + {Z1_A}) / 2 = {int_Z_A}")
print(f"Energy E_A = 1/2 * (Integral Z(x)) + 1 = 1/2 * {int_Z_A} + 1 = {E_A}")

print("\nCase B:")
print(f"Z(0) = T1(0) + T2(0) = {T1_0} + {T2_0_B} = {Z0_B}")
print(f"Z(1) = T1(sqrt(2)) + T2(1) = {T1_sqrt2} + {T2_1_B} = {Z1_B}")
print(f"Approximated Integral of Z(x) = (Z(0) + Z(1)) / 2 = ({Z0_B} + {Z1_B}) / 2 = {int_Z_B}")
print(f"Energy E_B = 1/2 * (Integral Z(x)) + 1 = 1/2 * {int_Z_B} + 1 = {E_B}")

print(f"\nThe minimum of the two energy values ({E_A} and {E_B}) is {min_E}.")
print(f"\nFinal Answer: E_total_min = 1/2 * (-1/2) + 1 = 3/4")
