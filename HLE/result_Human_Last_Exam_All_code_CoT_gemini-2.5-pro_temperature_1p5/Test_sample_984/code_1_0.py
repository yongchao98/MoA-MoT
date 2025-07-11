# Number of wells where paired heavy and light chains were successfully amplified
successful_wells = 4

# Total number of wells sorted with doublets (excluding the negative control)
total_sorted_wells = 95

# Calculate the efficiency as a decimal
efficiency_decimal = successful_wells / total_sorted_wells

# Convert the efficiency to a percentage
efficiency_percentage = efficiency_decimal * 100

# Print the final result in a clear format showing the equation
# The problem statement mentions low efficiency, let's quantify it.
print(f"The experiment yielded {successful_wells} successful wells out of {total_sorted_wells} sorted wells.")
print("The equation for the efficiency is:")
print(f"{successful_wells} / {total_sorted_wells} = {efficiency_decimal:.4f}")
print(f"This represents an efficiency of {efficiency_percentage:.2f}%.")

# Analysis of the low efficiency:
# The core challenge in this experiment is distinguishing true, antigen-specific cell interactions from random, non-specific cell clumps.
# - A & C: 40 minutes of co-incubation plus 30 minutes of staining on ice is a reasonable protocol. While not perfect, it is unlikely to be the primary cause for such a low success rate (around 4%).
# - D: The BD S8 is an imaging sorter specifically designed to improve doublet identification. Blaming the sorter for failing at its specialized task is less likely than a biological cause. Furthermore, if single B cells were sorted by mistake, you would still obtain heavy and light chains, they just wouldn't be specific to the tumor.
# - E: RNA competition from the larger tumor cell is a valid technical concern. However, it's a downstream molecular biology issue.
# - B: This addresses the fundamental biological principle at play. The frequency of B cells specific to any single antigen is inherently very low in a polyclonal population. Therefore, even if the sorter correctly identifies 100 physical doublets, the vast majority are expected to be random cell-cell stickiness rather than a stable, antigen-driven interaction. This means that most of the sorted "doublets" did not contain an antigen-specific B cell in the first place. This biological rarity is the most direct and compelling explanation for the very low number of positive hits.