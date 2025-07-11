import sympy

# Define symbols
r, gamma = sympy.symbols('r gamma')
xi = sympy.Function('xi')(r)
delta_P_elec = sympy.Function('Delta_P_elec')(r)

# The Young-Laplace equation relates the pressure difference (Delta_P) across an interface
# to the surface tension (gamma) and the mean curvature (H).
# Delta_P = 2 * gamma * H
# In three dimensions with two principal curvatures k1 and k2, Delta_P = gamma * (k1 + k2).

# For a surface of revolution z = xi(r) in cylindrical coordinates, the principal curvatures are:
# k1 = (d^2(xi)/dr^2) / (1 + (d(xi)/dr)^2)^(3/2)  (in the r-z plane)
# k2 = (d(xi)/dr) / (r * (1 + (d(xi)/dr)^2)^(1/2)) (azimuthal)

# Given the displacement is small, we use the small-slope approximation, where (d(xi)/dr)^2 << 1.
# This linearizes the curvature terms:
# k1_lin = d^2(xi)/dr^2
# k2_lin = (1/r) * d(xi)/dr

# The linearized Laplace pressure due to surface tension is:
# P_laplace = gamma * (k1_lin + k2_lin)
# P_laplace = gamma * (d^2(xi)/dr^2 + (1/r) * d(xi)/dr)

# The problem states that this pressure is balanced by the electrostatic pressure difference (Delta_P_elec).
# So, P_laplace = Delta_P_elec.

# The equation is: gamma * (d^2(xi)/dr^2 + (1/r) * d(xi)/dr) = Delta_P_elec
# Rearranging into the requested form: A(r)*xi'' + B(r)*xi' + C(r, xi) = 0
# gives: gamma*xi'' + (gamma/r)*xi' - Delta_P_elec = 0

# By comparing the terms, we can identify the coefficients A(r) and B(r).
A_r = gamma
B_r = gamma / r
C_r_xi = -delta_P_elec

# Display the final derived equation and the coefficients
print("The governing linear differential equation is:")
equation_str = f"{A_r} * d^2(xi)/dr^2 + ({B_r}) * d(xi)/dr - Delta_P_elec(r) = 0"
print(equation_str)
print("\nFrom this equation, we can identify the coefficients A(r) and B(r):")
print(f"A(r) = {A_r}")
print(f"B(r) = {B_r}")
