# Number of wells where paired heavy and light chains were successfully amplified
successful_wells = 4

# Total number of wells sorted with supposed doublets
total_sorted_wells = 95

# Calculate the experimental efficiency as a percentage
efficiency_percentage = (successful_wells / total_sorted_wells) * 100

# Print the calculation and the result
# The equation shows how the final efficiency percentage was derived.
print(f"The efficiency in obtaining paired heavy and light chains from sorted doublets was low.")
print(f"Calculation: ({successful_wells} successful wells / {total_sorted_wells} total sorted wells) * 100 = {efficiency_percentage:.2f}%")
print(f"This low success rate strongly suggests that most sorted doublets were not true antigen-specific pairs.")