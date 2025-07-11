import cmath

# BDF4 characteristic polynomial rho(xi)
def rho4(xi):
    # This is equivalent to the coefficient form:
    # (25/12)*xi**4 - 4*xi**3 + 3*xi**2 - (4/3)*xi + 1/4
    # It is derived from the sum formula which is more fundamental
    s = 0
    k = 4
    for j in range(1, k + 1):
        s += (1 / j) * (xi**(k - j)) * ((xi - 1)**j)
    return s

# Find the roots of N(xi) = xi^4 - (xi-1)^4 = 0
# The roots are (1+i)/2, (1-i)/2, 1/2
xi_cusp = (1 - 1j) / 2

# Value of rho at the cusp-defining xi
rho_val = rho4(xi_cusp)

# Value of xi^4 at that point
xi_pow4 = xi_cusp**4

# Location of the cusp on the stability boundary
z_cusp = rho_val / xi_pow4

# The angle alpha is determined by the tangent of the cusp's components
# alpha = arctan(|Im(z)| / |Re(z)|)
real_part = z_cusp.real
imag_part = z_cusp.imag

# The result is tan(alpha) = 8/3
print("The exact value of alpha is arctan(8/3).")
print("This can be written as:")
print("alpha = arctan({1}/{0})".format(abs(real_part), abs(imag_part)))
# We need to present the numbers in the final equation
y = abs(imag_part)
x = abs(real_part)
print("alpha = arctan( {} / {} )".format(y, x))
# Note that y/x = (2/3)/(1/4) = 8/3
numerator_of_ratio = 8
denominator_of_ratio = 3

print("So, the final value is arctan({}/{})".format(numerator_of_ratio, denominator_of_ratio))
