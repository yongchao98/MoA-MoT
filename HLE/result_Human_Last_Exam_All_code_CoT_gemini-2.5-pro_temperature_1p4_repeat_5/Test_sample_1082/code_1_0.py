# Define the sales price and environmental costs from the problem description.
sales_price = 50
water_cost = 10
energy_cost = 15

# Calculate the total environmental cost.
total_environmental_cost = water_cost + energy_cost

# Calculate the Sustainable Value Added (SVA).
sustainable_value_added = sales_price - total_environmental_cost

# Print the result in the format of a full equation to show the calculation.
print(f"Sustainable Value Added = ${sales_price} - (${water_cost} + ${energy_cost}) = ${sustainable_value_added}")