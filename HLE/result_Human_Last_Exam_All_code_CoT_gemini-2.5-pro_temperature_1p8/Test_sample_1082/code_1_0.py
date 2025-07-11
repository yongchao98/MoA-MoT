# Define the revenue and costs from the problem statement.
revenue = 50

# The phrase "...which cost 10 and 10 and 15 water and energy respectively"
# is interpreted as three separate costs for production.
cost_fruit = 10
cost_water = 10
cost_energy = 15

# Calculate the total cost of production.
total_cost = cost_fruit + cost_water + cost_energy

# Calculate the Sustainable Value Added (SVA).
sva = revenue - total_cost

# Print the equation step-by-step as requested.
print("To find the Sustainable Value Added (SVA), we subtract the total costs from the revenue.")
print(f"SVA = Revenue - (Cost of Fruits + Cost of Water + Cost of Energy)")
print(f"SVA = {revenue} - ({cost_fruit} + {cost_water} + {cost_energy})")
print(f"SVA = {revenue} - {total_cost}")
print(f"The final Sustainable Value Added is: {sva}")
