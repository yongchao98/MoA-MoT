# Define the vertices of the rectangle
rectangle_vertices = [(0, 0), (3, 0), (3, 2), (0, 2)]

# Calculate the area of the rectangle
rectangle_area = 3 * 2

# Define the vertices of the triangle where x < y
triangle_vertices = [(0, 0), (2, 2), (0, 2)]

# Calculate the area of the triangle
triangle_base = 2  # from (0, 2) to (2, 2)
triangle_height = 2  # from (0, 0) to (0, 2)
triangle_area = 0.5 * triangle_base * triangle_height

# Calculate the probability
probability = triangle_area / rectangle_area

print(probability)