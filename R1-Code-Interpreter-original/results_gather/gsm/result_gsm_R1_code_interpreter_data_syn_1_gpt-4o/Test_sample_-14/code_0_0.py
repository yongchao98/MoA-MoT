# Constants
average_weight_per_adult = 80  # kg
number_of_people = 9
maximum_load_capacity = 700  # kg

# Calculate total weight
total_weight = number_of_people * average_weight_per_adult

# Calculate excess weight
excess_weight = total_weight - maximum_load_capacity

# Output the result
print(excess_weight)