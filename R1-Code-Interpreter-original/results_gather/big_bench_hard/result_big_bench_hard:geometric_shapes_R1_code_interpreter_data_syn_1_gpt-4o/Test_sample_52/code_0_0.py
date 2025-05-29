# Define the path data
path_data = "M 29.40,69.39 L 51.16,67.08 L 65.98,58.97 L 86.05,37.54 L 83.22,21.48 M 83.22,21.48 L 41.67,26.57 L 29.40,69.39"

# Extract the coordinates from the path data
import re

# Use regex to find all coordinate pairs
coordinates = re.findall(r'(\d+\.\d+),(\d+\.\d+)', path_data)

# Convert to a set of unique points
unique_points = set(coordinates)

# Count the number of unique points
num_unique_points = len(unique_points)

# Output the number of unique points
print(num_unique_points)