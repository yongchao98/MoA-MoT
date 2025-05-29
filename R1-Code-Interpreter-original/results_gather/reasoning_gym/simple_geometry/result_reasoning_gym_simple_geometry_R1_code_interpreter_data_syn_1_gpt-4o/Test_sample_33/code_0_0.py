# Number of sides of the polygon
n = 6

# Sum of the interior angles of the polygon
sum_of_interior_angles = (n - 2) * 180

# Given angles
angles = [83.0, 83.0, 141.0, 66.0, 123.0]

# Calculate the sum of the given angles
sum_of_given_angles = sum(angles)

# Calculate the last angle
last_angle = sum_of_interior_angles - sum_of_given_angles

# Output the last angle
print(last_angle)