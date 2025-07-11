# The number of wells that successfully yielded paired heavy and light chains.
successful_wells = 4

# Total number of wells sorted into (96-well plate minus one negative control).
total_sorted_wells = 95

# Calculate the efficiency of the experiment.
efficiency_percentage = (successful_wells / total_sorted_wells) * 100

print(f"The experiment started with {total_sorted_wells} sorted wells and resulted in {successful_wells} successful outcomes.")
print(f"The success rate is calculated as: ({successful_wells} / {total_sorted_wells}) * 100 = {efficiency_percentage:.2f}%.")
print("\nThis low efficiency highlights a fundamental challenge in immunology.")
print("The most probable reason for this result is Choice B: True antigen-specific interactions are rare events, and many observed doublets are probably random associations rather than specific binding pairs.")
print("A cell sorter can identify that two cells are stuck together (a doublet), but it cannot distinguish between a true, high-affinity specific interaction and a random, non-specific clumping. Given the very low natural frequency of B cells specific for any single antigen, it is statistically likely that the vast majority of the sorted doublets were these random, non-specific pairings, which would not contain the correct B cell to yield a product in the subsequent PCR.")