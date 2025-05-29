# Given angles
angles = [14.0, 152.0, 144.0, 131.0, 24.0]

# Total sum of interior angles for a hexagon
total_sum = (6 - 2) * 180

# Sum of the given angles
sum_given_angles = sum(angles)

# Measure of the last angle
last_angle = total_sum - sum_given_angles

print(last_angle)