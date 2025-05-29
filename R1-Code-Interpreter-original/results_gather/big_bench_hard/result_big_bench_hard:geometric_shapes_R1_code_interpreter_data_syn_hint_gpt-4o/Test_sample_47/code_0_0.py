import math

def distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

# Points
points = [(41, 42), (37, 34), (41, 30), (45, 34), (41, 42)]

# Calculate distances
distances = [distance(points[i][0], points[i][1], points[i+1][0], points[i+1][1]) for i in range(len(points) - 1)]

print(distances)