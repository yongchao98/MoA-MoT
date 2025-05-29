# Given angles
angles = [79.0, 122.0, 127.0, 15.0, 101.0]

# Total sum of interior angles for a hexagon
total_sum = (6 - 2) * 180

# Sum of the given angles
sum_given_angles = sum(angles)

# Calculate the remaining angle
remaining_angle = total_sum - sum_given_angles

print(remaining_angle)