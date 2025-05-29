import math

# Given side lengths
a = 8
b = 4
c = 10

# Semi-perimeter
s = (a + b + c) / 2

# Area using Heron's formula
area_ABC = math.sqrt(s * (s - a) * (s - b) * (s - c))

# Area of triangle CDA
area_CDA = 25 * math.sqrt(3)

# Total area of quadrilateral ABCD
total_area = area_ABC + area_CDA

# Output the total area
print(total_area)