import collections

# Step 1: Establish the mappings based on the analysis.
# The keys are the plot letters (A-H) and values are the configuration numbers (1-8).
# This mapping is derived from matching the number of lobes/symmetry of the driver
# to the number of wiggles in the displacement plot, and then refining by shape.
mappings = {
    'A': 1,  # 3 complex wiggles -> Config 1 (3-fold complex symmetry)
    'B': 5,  # 6 small-amplitude wiggles -> Config 5 (6-fold smooth, smaller variation than 2)
    'C': 4,  # 4 sharp, high-amplitude wiggles -> Config 4 (4-fold sharp symmetry)
    'D': 8,  # 3 smooth wiggles -> Config 8 (3-fold simple symmetry)
    'E': 3,  # 1 major event/asymmetric curve -> Config 3 (asymmetric 2-lobe)
    'F': 7,  # 5 wiggles -> Config 7 (5-fold symmetry)
    'G': 6,  # 4 smooth wiggles -> Config 6 (4-fold smooth symmetry)
    'H': 2,  # 6 large-amplitude wiggles -> Config 2 (6-fold smooth, larger variation than 5)
}

# Sort the mappings alphabetically by plot letter to ensure correct order.
sorted_mappings = collections.OrderedDict(sorted(mappings.items()))

# Step 2: Print each individual pairing.
print("The determined pairings are:")
for plot, config in sorted_mappings.items():
    print(f"Plot {plot} corresponds to Configuration {config}")

# Step 3: Construct and print the final sequence.
final_sequence = "".join([str(config) for plot, config in sorted_mappings.items()])

print("\nThe final sequence of configuration numbers for plots A through H is:")
print(final_sequence)