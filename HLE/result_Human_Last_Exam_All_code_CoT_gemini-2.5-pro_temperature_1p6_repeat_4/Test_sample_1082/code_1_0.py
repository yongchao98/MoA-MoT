# This script calculates the Sustainable Value Added (SVA).

# 1. Define the variables based on the problem description.
# The final sale price of the product represents the economic value created.
sale_price = 50

# The costs for water and energy are based on the corporate eco-efficiency benchmark.
# These represent the environmental costs of the resources.
sustainable_water_cost = 10
sustainable_energy_cost = 15

# 2. Calculate the total environmental cost.
total_environmental_cost = sustainable_water_cost + sustainable_energy_cost

# 3. Calculate the Sustainable Value Added (SVA).
# SVA = Economic Value Created - Total Environmental Cost
sva = sale_price - total_environmental_cost

# 4. Print the final equation showing each number and the result.
# The format shows the value created minus the combined environmental costs.
print(f"Sustainable Value Added = {sale_price} - ({sustainable_water_cost} + {sustainable_energy_cost}) = {sva}")