# Define the financial variables from the problem description
revenue = 50  # The final product is sold for $50
cost_water = 10  # The cost of water is $10
cost_energy = 15  # The cost of energy is $15

# Calculate the total costs
total_costs = cost_water + cost_energy

# Calculate the Sustainable Value Added (SVA)
# SVA = Revenue - Total Costs
sva = revenue - total_costs

# Print the final equation with each number and the result
print(f"The Sustainable Value Added (SVA) is calculated as Revenue - (Water Cost + Energy Cost).")
print(f"SVA = {revenue} - ({cost_water} + {cost_energy})")
print(f"SVA = {revenue} - {total_costs}")
print(f"SVA = ${sva}")