# Define the financial variables based on the problem description.
# The sale price represents the revenue generated.
sale_price = 50

# The costs of water and energy represent the environmental costs.
cost_water = 10
cost_energy = 15

# Calculate the total environmental cost by summing the individual resource costs.
total_environmental_cost = cost_water + cost_energy

# Calculate the Sustainable Value Added (SVA) by subtracting the
# total environmental cost from the sale price.
sustainable_value_added = sale_price - total_environmental_cost

# Print the final equation showing all the components of the calculation.
print(f"Sustainable Value Added = Final Sale Price - (Cost of Water + Cost of Energy)")
print(f"Sustainable Value Added = ${sale_price} - (${cost_water} + ${cost_energy})")

# Print the final calculated value.
print(f"The sustainable value added for the product is: ${sustainable_value_added}")