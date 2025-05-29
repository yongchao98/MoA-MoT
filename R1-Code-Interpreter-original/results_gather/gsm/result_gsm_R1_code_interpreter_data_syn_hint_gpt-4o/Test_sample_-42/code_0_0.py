# Total value of equipment purchased
total_value = 4326221

# Percentage of faulty equipment
faulty_percentage = 40

# Calculate the value of faulty equipment
faulty_value = (faulty_percentage / 100) * total_value

# Calculate the value of functioning equipment
functioning_value = total_value - faulty_value

# Output the value of functioning equipment
print(functioning_value)