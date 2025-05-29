# Number of sides of the polygon
n = 6

# Sum of the interior angles of the polygon
total_sum_of_angles = (n - 2) * 180

# Given angles
angles = [56.0, 152.0, 65.0, 92.0, 85.0]

# Sum of the given angles
sum_of_given_angles = sum(angles)

# Measure of the last angle
last_angle = total_sum_of_angles - sum_of_given_angles

print(last_angle)