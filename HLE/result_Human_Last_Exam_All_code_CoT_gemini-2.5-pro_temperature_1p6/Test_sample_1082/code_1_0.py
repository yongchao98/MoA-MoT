# Define the financial and environmental cost values from the problem statement.
selling_price = 50
water_environmental_cost = 10
energy_environmental_cost = 15

# Calculate the total environmental cost by summing the individual costs.
total_environmental_cost = water_environmental_cost + energy_environmental_cost

# Calculate the Sustainable Value Added (SVA).
# SVA = Value Added (selling price) - Total Environmental Cost
sva = selling_price - total_environmental_cost

# Print the explanation and the final equation with the numbers plugged in.
print("To find the Sustainable Value Added (SVA), we subtract the total environmental costs from the product's selling price.")
print("\nFinal Equation:")
print(f"SVA = Selling Price - (Environmental Cost of Water + Environmental Cost of Energy)")
print(f"SVA = ${selling_price} - (${water_environmental_cost} + ${energy_environmental_cost}) = ${sva}")
