# Define the experimental parameters from the problem description
total_sorted_wells = 95
successful_wells_with_paired_chains = 4

# Calculate the experimental efficiency
# Efficiency is the ratio of successful outcomes to the total number of attempts, expressed as a percentage.
efficiency = (successful_wells_with_paired_chains / total_sorted_wells) * 100

# Print the final equation and the result
# The print statement explicitly shows each number used in the calculation.
print("The experiment aimed to isolate specific B cell:Tumor cell doublets from 95 sorted wells.")
print(f"Paired heavy and light chains were successfully amplified from {successful_wells_with_paired_chains} wells.")
print("The experimental efficiency is calculated as follows:")
print(f"({successful_wells_with_paired_chains} successful wells / {total_sorted_wells} total sorted wells) * 100 = {efficiency:.2f}%")