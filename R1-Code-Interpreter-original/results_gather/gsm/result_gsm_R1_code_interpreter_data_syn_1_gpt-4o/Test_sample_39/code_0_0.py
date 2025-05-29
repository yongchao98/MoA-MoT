# Constants
average_weight_per_adult = 80  # kg
number_of_adults = 9
maximum_load = 700  # kg

# Calculate total weight of the adults
total_weight = number_of_adults * average_weight_per_adult

# Calculate the excess weight
excess_weight = total_weight - maximum_load

# Output the result
print(excess_weight)