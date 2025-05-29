# Define the points from the SVG path
points = [
    (49.56, 30.61),
    (43.41, 40.44),
    (37.02, 36.44),
    (43.17, 26.61),
    (36.78, 22.62),
    (53.07, 18.45),
    (55.94, 34.61),
    (49.56, 30.61)  # Closing the path
]

# Count unique points
unique_points = set(points)
print(len(unique_points))