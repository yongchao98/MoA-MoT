# Known angles
angles = [127.0, 32.0, 42.0, 106.0, 18.0]

# Total sum of interior angles for a hexagon
total_sum = (6 - 2) * 180

# Sum of known angles
sum_known_angles = sum(angles)

# Remaining angle
remaining_angle = total_sum - sum_known_angles

print(remaining_angle)