# Number of wells where paired heavy and light chains were successfully amplified
successful_wells = 4

# Total number of wells sorted with presumed doublets
total_sorted_wells = 95

# Calculate the experimental efficiency as a percentage
efficiency = (successful_wells / total_sorted_wells) * 100

# Output the numbers used in the calculation and the result
print(f"Experimental Goal: Isolate antigen-specific Tumor:B-cell pairs.")
print(f"Successful wells with paired chains: {successful_wells}")
print(f"Total wells sorted with doublets: {total_sorted_wells}")
print(f"Efficiency Calculation: ({successful_wells} / {total_sorted_wells}) * 100")
print(f"Resulting Efficiency: {efficiency:.2f}%")