# Define the experimental parameters
total_sorted_wells = 95
successful_wells_with_paired_chains = 4

# Calculate the experimental efficiency
efficiency = (successful_wells_with_paired_chains / total_sorted_wells) * 100

# Output the equation and the result
print("Experimental Efficiency Calculation:")
print(f"({successful_wells_with_paired_chains} successful wells / {total_sorted_wells} total sorted wells) * 100 = {efficiency:.2f}%")
print("\n")
print("Analysis of Low Efficiency:")
print("The calculated efficiency is extremely low (~4.21%). This suggests a fundamental issue with the assumption that all sorted doublets represent specific interactions.")
print("The most likely reason for this low efficiency is:")
print("B. True antigen-specific interactions are rare events, and many observed doublets are probably random associations rather than specific binding pairs.")
print("\nExplanation:")
print("A cell sorter identifies doublets based on fluorescence signals (RFP+ and FITC+ appearing in the same event), not the underlying biological cause. While 95 doublets were sorted, it's highly probable that the vast majority were due to random, non-specific cell stickiness rather than a true, stable antigen-BCR interaction. Since the frequency of B cells specific to any single antigen is very low, it is expected that only a small fraction of the sorted B cells would be the correct ones.")
