# Define the financial values from the problem description.

# The selling price of the final product, which represents the Value Added.
selling_price = 50

# The cost of water used in the process.
water_cost = 10

# The cost of energy used in the process.
energy_cost = 15

# Calculate the total environmental impact cost.
environmental_impact_cost = water_cost + energy_cost

# Calculate the Sustainable Value Added (SVA).
sva = selling_price - environmental_impact_cost

# Print the formula and the final calculation step by step.
print("The formula for Sustainable Value Added (SVA) is: Value Added - Environmental Impact Cost")
print(f"In this case, SVA = Selling Price - (Water Cost + Energy Cost)")
print(f"SVA = ${selling_price} - (${water_cost} + ${energy_cost})")
print(f"SVA = ${selling_price} - ${environmental_impact_cost}")
print(f"The final Sustainable Value Added is: ${sva}")
