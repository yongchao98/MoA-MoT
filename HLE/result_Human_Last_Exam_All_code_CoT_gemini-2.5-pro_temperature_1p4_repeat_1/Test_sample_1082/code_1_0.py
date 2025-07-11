# Define the economic and environmental values from the problem description.

# The economic value is the selling price of the product.
economic_value = 50

# The environmental cost is interpreted as the total cost for water and energy used.
water_cost = 10
energy_cost = 15

# Calculate the total environmental cost.
total_environmental_cost = water_cost + energy_cost

# Calculate the Sustainable Value Added (SVA).
sva = economic_value - total_environmental_cost

# Print the final equation showing all the numbers involved in the calculation.
print(f"Sustainable Value Added = {economic_value} - ({water_cost} + {energy_cost}) = {sva}")