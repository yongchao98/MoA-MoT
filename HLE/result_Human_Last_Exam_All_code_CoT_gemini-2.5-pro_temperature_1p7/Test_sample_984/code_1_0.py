# Number of wells that were sorted into with presumed doublets.
# The total plate has 96 wells, but one was reserved for a negative control.
sorted_wells = 95

# Number of wells that successfully yielded paired heavy and light chains after PCR.
successful_wells = 4

# Calculate the efficiency of the experiment.
efficiency = (successful_wells / sorted_wells) * 100

print(f"Total number of sorted wells: {sorted_wells}")
print(f"Number of wells with successful paired chain amplification: {successful_wells}")
print(f"The final experimental efficiency was: {efficiency:.2f}%")
print("\nBased on this low efficiency, the most likely reason is that true antigen-specific interactions are rare events, and many of the sorted doublets were the result of random, non-specific cell association rather than true binding.")
