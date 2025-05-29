# Given perimeter of the square
perimeter_square = 108

# Calculate side length of the square
side_length_square = perimeter_square / 4

# Calculate the length of each segment after trisection
segment_length = side_length_square / 3

# Calculate the area of the square
area_square = side_length_square ** 2

# Calculate the area of one right-angled triangle
area_triangle = 0.5 * segment_length * segment_length

# Total area removed by four triangles
total_area_removed = 4 * area_triangle

# Calculate the area of the octagon
area_octagon = area_square - total_area_removed

print(area_octagon)