# Define the vertices of the hexagon
vertices = [(0, 0), (1, 3), (3, 3), (4, 0), (3, -3), (1, -3)]

# Initialize sums for the shoelace formula
sum1 = 0
sum2 = 0

# Calculate the sums
n = len(vertices)
for i in range(n):
    x1, y1 = vertices[i]
    x2, y2 = vertices[(i + 1) % n]  # Next vertex, wrapping around
    sum1 += x1 * y2
    sum2 += y1 * x2

# Calculate the area using the shoelace formula
area = 0.5 * abs(sum1 - sum2)

# Output the area
print(area)