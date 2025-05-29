# Given values
surface_area = 384
sum_of_edges = 112

# Calculate ab + bc + ca
ab_bc_ca = surface_area / 2

# Calculate a + b + c
a_b_c = sum_of_edges / 4

# Calculate (a + b + c)^2
a_b_c_squared = a_b_c ** 2

# Calculate a^2 + b^2 + c^2
a2_b2_c2 = a_b_c_squared - 2 * ab_bc_ca

# Calculate the radius r
import math
r = math.sqrt(a2_b2_c2) / 2

print(r)