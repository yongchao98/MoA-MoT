import re

# Path data
path_data = "M 34.08,19.65 L 43.45,30.98 L 38.19,35.33 L 28.82,24.00 L 26.80,25.67 L 24.63,13.16 L 36.11,17.98 M 36.11,17.98 L 34.08,19.65"

# Extract coordinates using regex
coordinates = re.findall(r'[\d.]+,[\d.]+', path_data)

# Convert to a set to get unique points
unique_points = set(coordinates)

# Count the number of unique points
num_unique_points = len(unique_points)

print(num_unique_points)