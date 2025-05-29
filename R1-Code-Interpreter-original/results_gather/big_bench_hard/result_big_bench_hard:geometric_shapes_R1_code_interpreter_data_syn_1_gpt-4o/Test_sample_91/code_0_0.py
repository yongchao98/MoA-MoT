# Define the path data
path_data = [
    (39.53, 40.66),
    (36.95, 20.64),
    (49.96, 18.96),
    (52.55, 38.98),
    (65.56, 37.30),
    (47.99, 58.90),
    (26.52, 42.34),
    (39.53, 40.66)
]

# Use a set to find unique points
unique_points = set(path_data)

# Count the number of unique points
num_unique_points = len(unique_points)

# Output the number of unique points
print(num_unique_points)