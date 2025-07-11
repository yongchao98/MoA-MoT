from fractions import Fraction

# Step 1: Start from the known root of the derived polynomial for x = 2*cos(theta)
# The relevant root for the tangent point is x = 2/5.
x_root = Fraction(2, 5)
# From x = 2*cos(theta), we find cos(theta)
c = x_root / 2

# Step 2: Calculate the real part (x_coord) of the tangent point z(theta).
# 12 * x_coord = 25 - 48*cos(theta) + 36*cos(2*theta) - 16*cos(3*theta) + 3*cos(4*theta)
# We use Chebyshev polynomials of the first kind (T_n) to express cos(n*theta) in terms of c = cos(theta)
cos_2theta = 2*c**2 - 1
cos_3theta = 4*c**3 - 3*c
cos_4theta = 8*c**4 - 8*c**2 + 1

# Calculate the numerator of 12 * x_coord
x_numerator_val = 25 - 48*c + 36*cos_2theta - 16*cos_3theta + 3*cos_4theta
# x_coord is this value divided by 12
x_coord = x_numerator_val / 12

# Step 3: Calculate the imaginary part (y_coord) of the tangent point z(theta).
# 12 * y_coord = 48*sin(theta) - 36*sin(2*theta) + 16*sin(3*theta) - 3*sin(4*theta)
# We factor out sin(theta) and use Chebyshev polynomials of the second kind (U_{n-1})
# where sin(n*theta)/sin(theta) = U_{n-1}(cos(theta)).
# 12 * y_coord / sin(theta) = 48*U_0(c) - 36*U_1(c) + 16*U_2(c) - 3*U_3(c)
# U_0(c)=1, U_1(c)=2c, U_2(c)=4c^2-1, U_3(c)=8c^3-4c
y_div_sin_numerator = 48 - 36*(2*c) + 16*(4*c**2 - 1) - 3*(8*c**3 - 4*c)
y_div_sin_coord = y_div_sin_numerator / 12

# Step 4: Calculate tan(alpha).
# The stability angle alpha is given by alpha = pi - arg(z).
# tan(alpha) = tan(pi - arg(z)) = -tan(arg(z)) = -y_coord / x_coord.
# tan_alpha = - (y_div_sin_coord * sin(theta)) / x_coord
# sin(theta) = sqrt(1 - c^2) = sqrt(1 - (1/5)^2) = sqrt(24/25) = 2*sqrt(6)/5.
# tan_alpha = -(y_div_sin_coord / x_coord) * (2*sqrt(6)/5)
rational_part = - y_div_sin_coord / x_coord
# The final result will be tan(alpha) = rational_part * 2 * sqrt(6) / 5
final_tan_alpha = rational_part * Fraction(2, 5)

A = final_tan_alpha.numerator
C = final_tan_alpha.denominator
B = 6 # The term under the square root

print("The A(alpha) stability angle alpha for BDF4 is given by alpha = arctan(A * sqrt(B) / C).")
print("The values are:")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print("\nSo the final equation is:")
print(f"alpha = arctan({A} * sqrt({B}) / {C})")