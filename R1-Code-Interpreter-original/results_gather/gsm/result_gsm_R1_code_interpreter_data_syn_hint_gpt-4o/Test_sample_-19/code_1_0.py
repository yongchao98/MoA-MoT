import math

# Given values
total_weight_kg = 2182736
weight_per_bag_kg = 50
cost_per_bag = 18

# Calculate the number of bags needed, rounding up to the nearest whole number
number_of_bags = math.ceil(total_weight_kg / weight_per_bag_kg)

# Calculate the total cost
total_cost = number_of_bags * cost_per_bag

# Print the total cost
print(total_cost)