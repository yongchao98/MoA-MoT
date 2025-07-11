import sympy

# Let's define the symbolic representation for the divergent part
# In the MS-bar scheme, divergences appear as poles in epsilon, where D = 4 - 2*epsilon.
# All one-loop counter-terms will be proportional to C = g**2 / (16 * sympy.pi**2 * epsilon)
# For simplicity in calculation, we will express all delta Z's in units of g**2 / (32 * pi**2 * epsilon)
# Let base_unit = g**2 / (32 * sympy.pi**2 * epsilon)
# This corresponds to C/2 in the line above.

# Step 1: Theoretical Calculation of Renormalization Constants

# Fermion self-energy Sigma_div = - (g**2 / (16 * pi**2 * epsilon)) * (p_slash/2 + M_x)
# Counter-term i(delta_Zx * p_slash - delta_m) cancels i*Sigma_div.
# So, delta_Zx * p_slash - delta_m = -Sigma_div = (g**2 / (16 * pi**2 * epsilon)) * (p_slash/2 + M_x)
#
# By comparing the coefficient of p_slash:
# delta_Zx = g**2 / (32 * pi**2 * epsilon)
# Expressed in our base unit, this is:
delta_Zx_coeff = 1

# By comparing the constant term (proportional to M_x):
# -delta_m = (g**2 / (16 * pi**2 * epsilon)) * M_x
# delta_m = - (g**2 / (16 * pi**2 * epsilon)) * M_x
#
# The mass renormalization constant is defined as delta_m = delta_Z_m_x * M_x
# So, delta_Z_m_x = delta_m / M_x = - g**2 / (16 * pi**2 * epsilon)
# Expressed in our base unit (which is half of this value), this is:
delta_Z_mx_coeff = -2

# The relation for the coupling renormalization constant is:
# Lambda_V_div / g = delta_Z_g + delta_Z_x + (1/2)*delta_Z_phi
# Given delta_Z_phi = 0, this simplifies to:
# Lambda_V_div / g = delta_Z_g + delta_Z_x
#
# The standard result for the vertex divergence is:
# Lambda_V_div / g = - g**2 / (16 * pi**2 * epsilon)
# So, delta_Z_g + delta_Z_x = - g**2 / (16 * pi**2 * epsilon)
#
# We can now solve for delta_Z_g:
# delta_Z_g = (-g**2 / (16 * pi**2 * epsilon)) - delta_Z_x
# delta_Z_g = (-g**2 / (16 * pi**2 * epsilon)) - (g**2 / (32 * pi**2 * epsilon))
# delta_Z_g = (-2 * g**2 / (32 * pi**2 * epsilon)) - (g**2 / (32 * pi**2 * epsilon))
# delta_Z_g = -3 * g**2 / (32 * pi**2 * epsilon)
# Expressed in our base unit, this is:
delta_Zg_coeff = -3

# Step 2: Calculate the ratio R

# R = delta_Z_x / (delta_Z_g + delta_Z_m_x)
# We can use the coefficients we found, as the common factor will cancel out.
R = delta_Zx_coeff / (delta_Zg_coeff + delta_Z_mx_coeff)

# Step 3: Print the results clearly

print("This script calculates the ratio R of one-loop renormalization constants in Yukawa theory.")
print("The calculation is performed in the MS-bar scheme with the condition delta_Z_phi = 0.")
print("\nLet C be the common factor g^2 / (32 * pi^2 * epsilon).\n")

print(f"1. Fermion field renormalization constant:")
print(f"   delta_Z_x = {delta_Zx_coeff} * C\n")

print(f"2. Fermion mass renormalization constant:")
print(f"   delta_Z_m_x = {delta_Z_mx_coeff} * C\n")

print(f"3. Yukawa coupling renormalization constant:")
print(f"   delta_Z_g = {delta_Zg_coeff} * C\n")

print("4. Now, we compute the ratio R = delta_Z_x / (delta_Z_g + delta_Z_m_x)")
print(f"   R = ({delta_Zx_coeff} * C) / (({delta_Zg_coeff} * C) + ({delta_Z_mx_coeff} * C))")
print(f"   R = {delta_Zx_coeff} / ({delta_Zg_coeff} + {delta_Z_mx_coeff})")
print(f"   R = {delta_Zx_coeff} / ({delta_Zg_coeff + delta_Z_mx_coeff})")
print(f"   R = {R}")
print("\nFinal answer for R:")
print(f"{R}")
print("\n<<<{}>>>".format(R))