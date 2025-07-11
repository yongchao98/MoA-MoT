# The number of wells where paired heavy and light chains were successfully amplified.
successful_wells = 4

# The total number of wells into which doublets were sorted.
total_sorted_wells = 95

# We calculate the efficiency by dividing the number of successful outcomes
# by the total number of attempts and multiplying by 100 to get a percentage.
efficiency_percentage = (successful_wells / total_sorted_wells) * 100

# We print the components of the calculation and the final result.
# The result is formatted to two decimal places for clarity.
print(f"Based on the experiment, the efficiency was calculated from {successful_wells} successful wells out of a total of {total_sorted_wells} sorted wells.")
print(f"Calculation: ({successful_wells} / {total_sorted_wells}) * 100 = {efficiency_percentage:.2f}%")
