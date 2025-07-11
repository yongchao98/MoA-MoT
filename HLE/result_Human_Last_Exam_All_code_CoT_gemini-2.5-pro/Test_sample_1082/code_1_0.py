# Step 1: Define the given variables based on the problem description.
selling_price = 50.0  # The final selling price of the product in dollars.
product_net_weight_g = 250.0  # The net weight of the product in grams.

# Step 2: Define the benchmark environmental costs.
# The costs are for an optimal/eco-efficient process.
# We interpret "cost 10and10and15 water and energy respectively" as a total
# cost of $10 for water and $15 for energy for a benchmark production batch.
benchmark_water_cost = 10.0
benchmark_energy_cost = 15.0

# Assume the benchmark cost corresponds to producing 1 kg (1000g) of product.
benchmark_production_weight_g = 1000.0

# Step 3: Calculate the total environmental cost for the benchmark unit (1 kg).
total_benchmark_cost = benchmark_water_cost + benchmark_energy_cost

# Step 4: Calculate the proportional environmental cost for our specific product (250g).
# This is the "sustainable" cost allowed for a product of this size.
product_environmental_cost = (product_net_weight_g / benchmark_production_weight_g) * total_benchmark_cost

# Step 5: Calculate the Sustainable Value Added (SVA).
# SVA = Value Created (Selling Price) - Environmental Cost Impact
sustainable_value_added = selling_price - product_environmental_cost

# Step 6: Print the final equation and the result.
# The output shows each number used in the final calculation step.
print("Calculating Sustainable Value Added (SVA):")
print(f"SVA = Selling Price - Proportional Environmental Cost")
print(f"SVA = ${selling_price:.2f} - ${product_environmental_cost:.2f}")
print(f"SVA = ${sustainable_value_added:.2f}")

# Final answer in the specified format
# print(f"<<<{sustainable_value_added}>>>")