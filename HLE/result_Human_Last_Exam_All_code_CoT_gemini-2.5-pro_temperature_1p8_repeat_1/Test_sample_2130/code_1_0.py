import math

# Plan:
# 1. The region traversed by the particles is a paraboloid of revolution.
# 2. We express its volume (V) and surface area (A) in terms of a single
#    dimensionless parameter, which we'll call 'y'.
#    The physical constraints of the problem require y > 1.
# 3. The ratio K = A^3 / V^2 can be written as a function of y:
#    K(y) = (16 * pi / 27) * (2*y^3 + 3*y^2 - 2)^3 / y^8
# 4. Analysis shows that K(y) is a monotonically increasing function for y > 1.
#    Therefore, the minimum value is the limit as y approaches its lower bound, 1.
# 5. We calculate this minimum value by substituting y = 1 into the expression.

print("The expression for the ratio to be minimized is K(y) = (16 * pi / 27) * (2*y^3 + 3*y^2 - 2)^3 / y^8.")
print("The minimum occurs at the boundary of the physical domain, which corresponds to y -> 1.")
print("\nWe evaluate the expression at y = 1:")
print("K_min = (16 * pi / 27) * (2*(1)^3 + 3*(1)^2 - 2)^3 / (1)^8")

# Define the numbers in the final equation
c1 = 16
c2 = 27
poly_term_coeff1 = 2
poly_term_val1 = 1
poly_term_pow1 = 3
poly_term_coeff2 = 3
poly_term_val2 = 1
poly_term_pow2 = 2
poly_term_const = 2
numerator_pow = 3
den_val = 1
den_pow = 8

print("\nDecomposition of the numbers in the equation:")
print(f"Constant factor: {c1} * pi / {c2}")
print(f"Polynomial term evaluated at y={poly_term_val1}:")
print(f"  {poly_term_coeff1}*({poly_term_val1}^{poly_term_pow1}) + {poly_term_coeff2}*({poly_term_val2}^{poly_term_pow2}) - {poly_term_const}")
print(f"The result of the polynomial is then raised to the power {numerator_pow}.")
print(f"Denominator term: {den_val}^{den_pow}")

# Calculate the value of the polynomial term
poly_val = poly_term_coeff1 * (poly_term_val1 ** poly_term_pow1) + poly_term_coeff2 * (poly_term_val2 ** poly_term_pow2) - poly_term_const

# Calculate the final minimum ratio
min_ratio = (c1 * math.pi / c2) * (poly_val ** numerator_pow) / (den_val ** den_pow)

print(f"\nCalculation steps:")
print(f"Value of the term in parentheses: 2*1 + 3*1 - 2 = {poly_val}")
print(f"This value is cubed: ({poly_val})^3 = {poly_val**numerator_pow}")
print(f"The ratio is: (16 * pi / 27) * {poly_val**numerator_pow} = 16 * pi")

print("\nFinal numerical answer:")
print(min_ratio)
<<<50.26548245743669>>>