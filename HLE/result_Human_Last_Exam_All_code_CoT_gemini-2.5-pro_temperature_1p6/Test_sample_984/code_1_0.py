# Define the experimental parameters based on the problem description.
# The total number of wells into which doublets were sorted.
total_sorted_wells = 95

# The number of wells that successfully yielded paired heavy and light chains after PCR.
successful_wells = 4

# Calculate the efficiency of the sorting and amplification process.
# The formula is (successful outcomes / total attempts) * 100.
efficiency_percentage = (successful_wells / total_sorted_wells) * 100

# Output a sentence that clearly states the numbers involved and the final calculated efficiency.
# This quantifies the "low efficiency" described in the problem.
print(
    f"Out of {total_sorted_wells} sorted wells, only {successful_wells} "
    f"yielded the desired paired chains, which represents a "
    f"success rate of {efficiency_percentage:.2f}%."
)