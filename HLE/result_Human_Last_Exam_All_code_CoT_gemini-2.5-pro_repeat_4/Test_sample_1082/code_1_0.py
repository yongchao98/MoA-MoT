# Step 1: Define the economic value and environmental costs from the problem description.
selling_price = 50
cost_of_water = 10
cost_of_energy = 15

# Step 2: Calculate the total environmental cost.
total_environmental_cost = cost_of_water + cost_of_energy

# Step 3: Calculate the Sustainable Value Added (SVA).
# SVA is the value created (selling price) minus the environmental costs.
sustainable_value_added = selling_price - total_environmental_cost

# Step 4: Print the explanation, the final equation with all its components, and the result.
print("To calculate the Sustainable Value Added (SVA), we subtract the total environmental costs from the product's selling price.")
print("The final equation is:")
print(f"{selling_price} (Selling Price) - ({cost_of_water} (Water Cost) + {cost_of_energy} (Energy Cost)) = {sustainable_value_added}")
print(f"\nThe Sustainable Value Added for the product is ${sustainable_value_added}.")