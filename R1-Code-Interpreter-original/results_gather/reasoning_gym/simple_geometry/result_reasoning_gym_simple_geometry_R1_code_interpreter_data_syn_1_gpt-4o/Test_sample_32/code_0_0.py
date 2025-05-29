# Given angles
angles = [41.0, 166.0, 113.0, 147.0, 14.0]

# Total sum of interior angles for a hexagon
total_sum = (6 - 2) * 180

# Sum of the given angles
sum_given_angles = sum(angles)

# Calculate the remaining angle
remaining_angle = total_sum - sum_given_angles

print(remaining_angle)