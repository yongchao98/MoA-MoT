import sympy as sp

# Define symbols
R, H = sp.symbols('R H', positive=True)

# Heights of each section
h = H / 4

# Radii of the slices
r1 = R
r2 = 3*R/4
r3 = R/2

# Volume of the largest piece (bottom frustum)
V_largest = (1/3) * sp.pi * h * (r1**2 + r1*r2 + r2**2)

# Volume of the second-largest piece
V_second_largest = (1/3) * sp.pi * h * (r2**2 + r2*r3 + r3**2)

# Ratio of the volumes
ratio = V_second_largest / V_largest

# Simplify the ratio
simplified_ratio = sp.simplify(ratio)

print(simplified_ratio)