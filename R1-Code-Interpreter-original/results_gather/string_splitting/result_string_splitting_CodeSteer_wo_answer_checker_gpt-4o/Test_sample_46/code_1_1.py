# Initial counts of machines
machines_A = 2
machines_B = 2
machines_C = 2

# Calculate maximum parts X and Y from machines A, B, and C
max_parts_X_from_A = machines_A * 2
max_parts_Y_from_A = machines_A

max_parts_X_from_B = (machines_B // 2)
max_parts_Y_from_C = (machines_C // 2)

# Total possible parts X and Y
total_parts_X = max_parts_X_from_A + max_parts_X_from_B
total_parts_Y = max_parts_Y_from_A + max_parts_Y_from_C

# Calculate maximum parts Z that can be formed
max_parts_Z = min(total_parts_X, total_parts_Y)

# Remaining parts X and Y after forming parts Z
remaining_parts_X = total_parts_X - max_parts_Z
remaining_parts_Y = total_parts_Y - max_parts_Z

# Calculate remaining machines
remaining_machines_A = machines_A - (max_parts_X_from_A // 2)
remaining_machines_B = machines_B % 2
remaining_machines_C = machines_C % 2

# Output the final counts
result = [
    str(remaining_machines_A),
    str(remaining_machines_B),
    str(remaining_machines_C),
    str(remaining_parts_X),
    str(remaining_parts_Y),
    str(max_parts_Z)
]
print(result)