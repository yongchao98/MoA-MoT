import math

# Given values
total_weight = 2182736  # in kg
weight_per_bag = 50  # in kg
cost_per_bag = 18  # in dollars

# Calculate the number of bags, rounding up to ensure all coal is covered
number_of_bags = math.ceil(total_weight / weight_per_bag)

# Calculate the total cost
total_cost = number_of_bags * cost_per_bag

# Print the total cost
print(total_cost)