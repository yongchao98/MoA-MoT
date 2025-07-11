# Define the financial values from the problem statement.
revenue = 50
water_cost = 10
energy_cost = 15

# Calculate the total environmental cost by summing the cost of water and energy.
total_environmental_cost = water_cost + energy_cost

# Calculate the Sustainable Value Added (SVA).
sva = revenue - total_environmental_cost

# Print the final equation with all the numbers, as requested.
print(f"Sustainable Value Added (SVA) = ${revenue} (Revenue) - (${water_cost} (Water Cost) + ${energy_cost} (Energy Cost))")

# Print the final result.
print(f"SVA = ${sva}")