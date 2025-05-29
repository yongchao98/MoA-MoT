import sympy as sp

# Define symbols
R, H = sp.symbols('R H')

# Volume of the largest piece (bottom frustum)
V_largest = (1/3) * sp.pi * (H/4) * ((3*R/4)**2 + (3*R/4)*R + R**2)

# Volume of the second-largest piece (third frustum)
V_second_largest = (1/3) * sp.pi * (H/4) * ((R/2)**2 + (R/2)*(3*R/4) + (3*R/4)**2)

# Calculate the ratio
ratio = V_second_largest / V_largest

# Simplify the ratio
simplified_ratio = sp.simplify(ratio)

# Express the ratio as a common fraction
common_fraction = sp.Rational(simplified_ratio).limit_denominator()

print(common_fraction)