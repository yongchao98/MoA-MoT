# The number of wells from which paired heavy and light chains were successfully recovered.
successful_wells = 4

# The total number of wells that received a sorted doublet.
sorted_wells = 95

# The efficiency of the experiment is the ratio of successful outcomes to the total number of attempts.
efficiency = successful_wells / sorted_wells

# To fulfill the request, we print the components of the equation and the final result.
print("The efficiency of recovering antigen-specific paired chains is calculated as:")
print(f"{successful_wells} (successful wells) / {sorted_wells} (total sorted wells) = {efficiency}")
print(f"This is an efficiency of approximately {efficiency:.1%}.")