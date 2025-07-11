# Define the financial variables based on the problem description.
# The interpretation assumes the costs provided are totals for this production process.
product_revenue = 50
raw_material_cost = 10 # Total cost for fresh mangoes and coconuts.
water_cost = 10 # Total cost for water consumption.
energy_cost = 15 # Total cost for energy consumption.

# Calculate the economic value added before considering environmental costs.
value_added = product_revenue - raw_material_cost

# Calculate the total environmental cost.
total_environmental_cost = water_cost + energy_cost

# Calculate the final Sustainable Value Added (SVA).
sustainable_value_added = value_added - total_environmental_cost

# Output the step-by-step calculation and the final equation with all numbers.
print("Sustainable Value Added Calculation:")
print("Formula: SVA = (Product Revenue - Raw Material Cost) - (Water Cost + Energy Cost)")
print("\nStep 1: Calculate Economic Value Added")
print(f"Value Added = ${product_revenue} (Revenue) - ${raw_material_cost} (Raw Materials) = ${value_added}")

print("\nStep 2: Calculate Total Environmental Cost")
print(f"Environmental Cost = ${water_cost} (Water) + ${energy_cost} (Energy) = ${total_environmental_cost}")

print("\nStep 3: Calculate Sustainable Value Added (SVA)")
print("Final Equation:")
print(f"SVA = ${value_added} (Value Added) - ${total_environmental_cost} (Environmental Cost)")
print(f"SVA = (${product_revenue} - ${raw_material_cost}) - (${water_cost} + ${energy_cost}) = ${sustainable_value_added}")

print(f"\nThe final Sustainable Value Added is: ${sustainable_value_added}")