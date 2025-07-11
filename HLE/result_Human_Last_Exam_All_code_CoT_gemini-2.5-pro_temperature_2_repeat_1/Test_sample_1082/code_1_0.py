# Define the financial values provided in the problem description
revenue = 50
water_cost = 10
energy_cost = 15

# Calculate the total costs of the resources used
total_cost = water_cost + energy_cost

# Calculate the Sustainable Value Added (SVA).
# Since the process is described as eco-efficient ("optimal use"),
# the value added (Revenue - Costs) is the Sustainable Value Added.
sva = revenue - total_cost

# Print the final equation with all the numbers, as requested.
print(f"The Sustainable Value Added (SVA) is calculated as Revenue minus the sum of resource costs.")
print(f"SVA = ${revenue} - (${water_cost} + ${energy_cost}) = ${sva}")