# Number of wells where paired heavy and light chains were successfully amplified
successful_wells = 4

# Total number of wells that were sorted with presumed doublets
total_sorted_wells = 95

# Calculate the experimental efficiency as a percentage
efficiency = (successful_wells / total_sorted_wells) * 100

# Print the equation and the resulting efficiency
print("The efficiency of recovering specific paired chains is calculated as follows:")
print(f"Equation: ({successful_wells} successful wells / {total_sorted_wells} total sorted wells) * 100")
print(f"Result: {efficiency:.2f}%")
print("\nThis low efficiency strongly suggests that the vast majority of sorted 'doublets' were not true antigen-specific pairs.")
