import sympy

# Define symbols for the parameters.
g, epsilon = sympy.symbols('g epsilon')
# In MS-bar scheme, counter-terms are proportional to 1/epsilon.
# The overall constant C captures the common factors.
C = g**2 / (16 * sympy.pi**2 * epsilon)

# From the fermion self-energy calculation, we have the following results for the counter-terms:
# 1. Fermion wave function renormalization counter-term
delta_Zx = C / 2

# 2. A combination of fermion mass and wave function counter-terms
sum_delta_Z_xm = -C

# From the Ward Identity for Yukawa coupling, we assume delta_Zg = delta_Zx.
delta_Zg = delta_Zx

# From the sum_delta_Z_xm expression, we can find delta_Zm_x.
# sum_delta_Z_xm = delta_Zx + delta_Zm_x
# so, delta_Zm_x = sum_delta_Z_xm - delta_Zx
delta_Zm_x = sum_delta_Z_xm - delta_Zx

# Now, we calculate the ratio R = delta_Zx / (delta_Zg + delta_Zm_x)
# Substitute the expression for delta_Zg and delta_Zm_x
denominator = delta_Zg + delta_Zm_x

R = delta_Zx / denominator

# Print the step-by-step derivation
print("Calculation of the ratio R")
print("==========================")
print(f"1. The fermion field counter-term is: δZ_x = {delta_Zx}")
print(f"2. The sum of δZ_x and δZ_m_x is: δZ_x + δZ_m_x = {sum_delta_Z_xm}")
print(f"3. Using the Ward Identity (δZ_g = δZ_x), the Yukawa coupling counter-term is: δZ_g = {delta_Zg}")
print(f"4. From steps 1 and 2, we can find δZ_m_x = (δZ_x + δZ_m_x) - δZ_x = {sum_delta_Z_xm} - ({delta_Zx}) = {delta_Zm_x}")
print(f"5. The denominator of R is δZ_g + δZ_m_x = ({delta_Zg}) + ({delta_Zm_x}) = {denominator}")
print(f"6. The ratio R is δZ_x / (δZ_g + δZ_m_x).")
print(f"   R = ({delta_Zx}) / ({denominator})")
print(f"\nFinal Equation: R = (g^2 / (32 * π^2 * ε)) / -(g^2 / (16 * π^2 * ε)) = {R}")