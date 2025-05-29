import math

# Side length of the octagon
s = 1  # Assume side length is 1 for simplicity

# Area of the regular octagon
A_octagon = 2 * (1 + math.sqrt(2)) * s**2

# Area of one isosceles triangle (like AOB)
A_triangle = A_octagon / 8

# Area of ABMO
# ABMO = AOB + BMO
# BMO is half of BCO, which is another isosceles triangle
A_ABMO = A_triangle + (A_triangle / 2)

# Area of EDCMO
# EDCMO = EDCM + MCO
# MCO is half of BCO
# EDCM is a quadrilateral, which is 2 triangles (EDO + DCO)
A_EDCMO = 2 * A_triangle + (A_triangle / 2)

# Ratio of the areas
ratio = A_ABMO / A_EDCMO

print(ratio)