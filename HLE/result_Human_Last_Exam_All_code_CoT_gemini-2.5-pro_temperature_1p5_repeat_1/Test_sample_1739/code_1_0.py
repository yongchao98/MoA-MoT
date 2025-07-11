# The nonlinear correction to the linear oscillation frequency is given by a polynomial
# in terms of the polytropic index, gamma.
# The polynomial part is derived from the expression: 2*omega_0^4 - 3*omega_0^2 - 6,
# where omega_0^2 = 3*gamma.

# Let the polynomial be A*gamma^2 + B*gamma + C.
# We calculate the coefficients A, B, and C.

# Coefficient of the gamma^2 term (from 2*omega_0^4 = 2*(3*gamma)^2)
coeff_A = 2 * (3**2)

# Coefficient of the gamma term (from -3*omega_0^2 = -3*(3*gamma))
coeff_B = -3 * 3

# The constant term
coeff_C = -6

print("The polynomial factor in the nonlinear frequency correction is of the form: A*gamma^2 + B*gamma + C")
print(f"The equation with the calculated numbers is: {coeff_A}*gamma^2 + ({coeff_B})*gamma + ({coeff_C})")
print("\nThe terms of the correction polynomial are:")
print(f"1st term: {coeff_A}*gamma^2")
print(f"2nd term: {coeff_B}*gamma")
print(f"3rd term: {coeff_C}")
