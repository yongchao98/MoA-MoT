# Define the numbers provided in the experimental description
successful_wells = 4
total_sorted_wells = 95

# Calculate the experimental efficiency as a percentage
efficiency = (successful_wells / total_sorted_wells) * 100

# Output the components of the calculation and the final result
print("This script calculates the experimental efficiency.")
print(f"Number of wells with paired heavy and light chains: {successful_wells}")
print(f"Total number of wells sorted with doublets: {total_sorted_wells}")
print("The efficiency calculation is:")
print(f"({successful_wells} / {total_sorted_wells}) * 100 = {efficiency:.2f}%")