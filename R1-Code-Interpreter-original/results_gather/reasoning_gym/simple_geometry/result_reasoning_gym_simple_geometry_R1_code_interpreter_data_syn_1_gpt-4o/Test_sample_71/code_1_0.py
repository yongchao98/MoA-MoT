# Given angles
angles = [65.0, 29.0, 164.0, 41.0, 107.0]

# Total sum of interior angles for a hexagon
total_sum = (6 - 2) * 180

# Sum of the given angles
sum_of_given_angles = sum(angles)

# Calculate the last angle
last_angle = total_sum - sum_of_given_angles

print("Sum of given angles:", sum_of_given_angles)
print("Last angle:", last_angle)