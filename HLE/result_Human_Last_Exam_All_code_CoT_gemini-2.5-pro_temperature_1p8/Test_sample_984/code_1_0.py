# The user wants to understand the low efficiency of their experiment.
# The calculation for the efficiency is based on the numbers provided.

# Number of wells with paired heavy and light chains successfully amplified.
successful_wells = 4

# Total number of sorted wells (excluding the negative control).
total_sorted_wells = 95

# Calculate the efficiency.
efficiency = (successful_wells / total_sorted_wells) * 100

# The core of the problem is to choose the best explanation from the provided options.
# Based on the step-by-step thinking process, the most fundamental reason for the
# low efficiency is the biological rarity of the specific interaction. Most of the
# sorted "doublets" are not the real target, but are random, weak associations
# that do not survive the sorting process, leading to PCR failure.

# This block prints out the analysis in text form.
# The final answer will be selected from the multiple-choice options.

print("Analysis of Experimental Efficiency:")
print(f"Number of successful wells: {successful_wells}")
print(f"Total sorted wells: {total_sorted_wells}")
print(f"Experimental efficiency: {successful_wells}/{total_sorted_wells} = {efficiency:.2f}%")
print("\nConclusion:")
print("The extremely low efficiency (~4.2%) is best explained by the fundamental biological principle that true antigen-specific interactions are rare events.")
print("Most of the doublets detected by the sorter are likely due to random, non-specific cell associations.")
print("These random associations are much weaker than true specific binding and are unlikely to survive the physical shear forces of the sorting process.")
print("As a result, most wells do not receive the intact B cell:tumor cell pair required for successful downstream amplification.")
print("Therefore, the low efficiency in obtaining paired chains is a direct consequence of the rarity of stable, specific doublets.")
