# Dimensions of the rectangle
width = 2008
height = 2009

# Area of the rectangle
area_rectangle = width * height

# Base and height of the triangle
base_triangle = width
height_triangle = width / 2  # 1004

# Area of the triangle
area_triangle = 0.5 * base_triangle * height_triangle

# Probability
probability = area_triangle / area_rectangle

# Output the probability as a common fraction
from fractions import Fraction
probability_fraction = Fraction(probability).limit_denominator()

print(f"<<<{probability_fraction}>>>")