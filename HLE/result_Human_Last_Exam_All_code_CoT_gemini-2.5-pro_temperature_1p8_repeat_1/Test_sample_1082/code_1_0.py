# Step 1: Define the financial values from the problem description.
revenue = 50  # The product is sold for $50.
water_cost = 10  # The cost attributed to water usage.
energy_cost = 15  # The cost attributed to energy usage.

# Step 2: Calculate the total environmental cost.
total_environmental_cost = water_cost + energy_cost

# Step 3: Calculate the Sustainable Value Added (SVA).
sustainable_value_added = revenue - total_environmental_cost

# Step 4: Print the final equation and the result.
# The calculation shows the revenue minus the sum of the environmental costs.
print("Sustainable Value Added Calculation:")
print(f"{revenue} (Revenue) - ({water_cost} (Water Cost) + {energy_cost} (Energy Cost)) = {sustainable_value_added}")