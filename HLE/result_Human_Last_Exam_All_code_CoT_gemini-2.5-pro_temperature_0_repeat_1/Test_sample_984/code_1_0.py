# Number of wells where paired heavy and light chains were successfully amplified
successful_wells = 4

# Total number of wells into which doublets were sorted
total_sorted_wells = 95

# Calculate the efficiency as a percentage
efficiency = (successful_wells / total_sorted_wells) * 100

# Print the final equation and the result
print(f"The experimental efficiency is calculated as ({successful_wells} / {total_sorted_wells}) * 100.")
print(f"Efficiency: {efficiency:.2f}%")