# Define the dimensions of the rectangle
rectangle_width = 2008
rectangle_height = 2009

# Calculate the area of the rectangle
rectangle_area = rectangle_width * rectangle_height

# Define the base and height of the triangle
triangle_base = 2008
triangle_height = 1004

# Calculate the area of the triangle
triangle_area = 0.5 * triangle_base * triangle_height

# Calculate the probability
probability = triangle_area / rectangle_area

# Print the probability as a common fraction
from fractions import Fraction
print(Fraction(probability).limit_denominator())