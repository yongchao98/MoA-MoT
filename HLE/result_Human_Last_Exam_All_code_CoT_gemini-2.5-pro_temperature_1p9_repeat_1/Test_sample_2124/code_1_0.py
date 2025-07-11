# Let's define the values of the counter-term coefficients
# in units of the common factor C = g^2 / (32 * pi^2 * epsilon).
# According to our derivation:
# delta_Zx = -1 * C
# delta_Zmx = 2 * C
# delta_Zg = -3 * C
dZ_x_coeff = -1
dZ_mx_coeff = 2
dZ_g_coeff = -3

# Now, we calculate the ratio R.
# The common factor C will cancel out.
R_numerator = dZ_x_coeff
R_denominator = dZ_g_coeff + dZ_mx_coeff

# To avoid division by zero, though in this case the denominator is not zero.
if R_denominator == 0:
    print("The denominator is zero, the ratio R is undefined.")
else:
    R = R_numerator / R_denominator
    print(f"The numerator is the coefficient of delta_Zx: {R_numerator}")
    print(f"The denominator is the sum of coefficients of delta_Zg and delta_Zmx: {dZ_g_coeff} + {dZ_mx_coeff} = {R_denominator}")
    print(f"The ratio R is therefore {R_numerator} / {R_denominator} = {R}")
