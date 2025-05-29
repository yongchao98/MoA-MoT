# Given angles
angles = [138.0, 38.0, 165.0, 158.0, 72.0]

# Total sum of interior angles for a hexagon
total_sum = (6 - 2) * 180

# Sum of the given angles
sum_of_given_angles = sum(angles)

# Calculate the last angle
last_angle = total_sum - sum_of_given_angles

# Output the last angle
print(last_angle)