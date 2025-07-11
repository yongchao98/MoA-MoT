# The problem is to identify the most likely reason for a low experimental efficiency.
# Let's first calculate the exact efficiency to quantify how low it is.

# Total number of wells sorted with presumed doublets.
total_sorted_wells = 95

# Number of wells that successfully yielded paired heavy and light chains.
successful_wells = 4

# Calculate the efficiency percentage.
efficiency = (successful_wells / total_sorted_wells) * 100

# Print the numbers used in the calculation and the final result.
print(f"From a total of {total_sorted_wells} sorted wells, paired heavy and light chains were recovered from {successful_wells} wells.")
print(f"The final efficiency is calculated as: {successful_wells} / {total_sorted_wells} * 100 = {efficiency:.2f}%")

print("\nAnalysis:")
print("This extremely low efficiency of ~4.2% suggests a fundamental issue with the experimental assumptions.")
print("The most likely reason is that while 95 'doublet' events were sorted based on fluorescence, the vast majority were not true, stable antigen-specific interactions.")
print("True antigen-specific B cells are rare, and many cells will stick together randomly. The sorter cannot distinguish between a specific biological interaction and a random, non-specific aggregate. Therefore, it is expected that most sorted doublets were random pairings, leading to a low success rate in recovering chains from genuinely interacting B cells.")