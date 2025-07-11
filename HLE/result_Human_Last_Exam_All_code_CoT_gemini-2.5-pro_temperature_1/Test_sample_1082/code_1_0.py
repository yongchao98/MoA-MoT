import sys

# --- Input Values from the problem statement ---
product_net_weight_g = 250.0
# The ratio of mangoes to coconuts (1:2) is descriptive information about the product's composition.
# It is not needed for the SVA calculation itself.

# --- Benchmark Production Batch Information ---
total_revenue = 50.0  # in dollars
water_cost = 10.0  # in dollars
energy_cost = 15.0  # in dollars
water_used_liters = 100.0

# This is interpreted as the yield of the dehydration process, where
# final dry weight = 0.20 * initial fresh weight.
yield_factor = 0.20

# Step 1: Calculate the total economic value added for the batch
print("Step 1: Calculate the total economic value added for the production batch.")
total_resource_cost = water_cost + energy_cost
total_value_added = total_revenue - total_resource_cost
print(f"Total Revenue: ${total_revenue:.2f}")
print(f"Total Resource Cost (Water + Energy): ${water_cost:.2f} + ${energy_cost:.2f} = ${total_resource_cost:.2f}")
print(f"Total Value Added for the batch: ${total_revenue:.2f} - ${total_resource_cost:.2f} = ${total_value_added:.2f}\n")

# Step 2: Determine the total weight of the production batch
print("Step 2: Determine the total weight of dry fruit produced in the batch.")
# We assume the 100 L of water used is the water removed from the fruit during dehydration.
# 1 Liter of water has a mass of approximately 1 kg.
water_removed_kg = water_used_liters

# The mass of water removed is the difference between the initial fresh weight and the final dry weight.
# Water Removed = Fresh Weight - Dry Weight
# Fresh Weight = Dry Weight / Yield Factor
# Water Removed = (Dry Weight / Yield Factor) - Dry Weight
# Water Removed = Dry Weight * (1 / Yield Factor - 1)
# So, Dry Weight = Water Removed / (1 / Yield Factor - 1)
total_dry_weight_kg = water_removed_kg / ((1 / yield_factor) - 1)
total_dry_weight_g = total_dry_weight_kg * 1000

print(f"The calculation for total batch weight is based on the amount of water removed ({water_removed_kg:.0f} kg) and the fruit yield ({yield_factor*100:.0f}%).")
print(f"Total Batch Weight (Dry Fruit) = {water_removed_kg:.0f} kg / ((1 / {yield_factor}) - 1) = {total_dry_weight_kg:.0f} kg or {total_dry_weight_g:.0f} g\n")


# Step 3: Calculate the Sustainable Value Added (SVA) for the 250g product
print("Step 3: Calculate the Sustainable Value Added (SVA) for the specific 250g product.")
# SVA is the product's proportional share of the total value added.
sva_for_product = (product_net_weight_g / total_dry_weight_g) * total_value_added
print(f"The SVA is the product's share of the total value added.")
print(f"SVA = (Product Weight / Total Batch Weight) * Total Value Added")
# Printing the final equation with all numbers
equation = f"SVA = ({product_net_weight_g:.0f}g / {total_dry_weight_g:.0f}g) * (${total_revenue:.0f} - (${water_cost:.0f} + ${energy_cost:.0f}))"
print("Final Equation:")
print(equation)
print(f"SVA = {product_net_weight_g/total_dry_weight_g} * ${total_value_added:.0f} = ${sva_for_product:.2f}")

# The final answer is the numerical value of the SVA
sys.stdout.write("\n")
sys.stdout.write(f'<<<{sva_for_product:.2f}>>>')