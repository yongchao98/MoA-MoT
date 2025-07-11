# Plan:
# 1. Define the number of successful wells.
# 2. Define the total number of wells where sorting occurred.
# 3. Calculate the efficiency percentage.
# 4. Print the result in a clear, explanatory sentence, showing the original numbers used in the calculation.

# The problem states that out of 95 sorted wells, 4 yielded the correct paired chains.
# Let's calculate the experimental efficiency.

successful_wells = 4
total_sorted_wells = 95

# Calculate the efficiency as a percentage
efficiency = (successful_wells / total_sorted_wells) * 100

# Print the final result in a sentence.
# The final answer to the multiple-choice question is B.
# The code below serves to quantify the "low efficiency" mentioned.

print(f"The experimental efficiency was calculated by dividing the number of successful wells ({successful_wells}) by the total number of sorted wells ({total_sorted_wells}) and multiplying by 100.")
print(f"Calculation: ({successful_wells} / {total_sorted_wells}) * 100 = {efficiency:.2f}%")
print("\nThis extremely low efficiency is best explained by the biological rarity of the event.")
print("The most likely reason is that true antigen-specific interactions are rare, and most of the sorted doublets were the result of random, non-specific cell aggregation rather than specific binding.")
