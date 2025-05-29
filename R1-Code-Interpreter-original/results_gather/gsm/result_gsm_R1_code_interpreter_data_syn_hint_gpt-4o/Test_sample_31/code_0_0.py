# Total number of trucks
total_trucks = 1088322

# Tons of fertiliser per truck
tons_per_truck = 1088322

# Calculate the number of trucks that failed
failed_trucks = total_trucks // 4

# Calculate the number of successful trucks
successful_trucks = total_trucks - failed_trucks

# Calculate the total tons of fertiliser delivered
total_tons_delivered = successful_trucks * tons_per_truck

print(total_tons_delivered)