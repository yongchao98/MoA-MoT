import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Step 1: Define the input parameters based on the problem statement ---

# Product specific information
product_weight_g = 250.0  # Net weight of the final product in grams

# Batch process information
# We interpret "100 liters of water" as the amount of water removed during dehydration.
# Assuming 1 liter of water has a mass of 1 kg.
water_removed_kg = 100.0
energy_used_kWh = 40.0  # Energy consumption for the batch
moisture_content_fresh = 0.20  # 20% moisture in fresh fruits

# Batch financial information
total_revenue = 50.0  # Revenue from selling the entire batch
# Interpreting the ambiguous phrase "cost 10 and 10 and 15 water and energy respectively"
# as costs for Fruit (raw material), Water, and Energy for the batch.
cost_fruit = 10.0
cost_water = 10.0
cost_energy = 15.0

# --- Step 2: Calculate the total weight of the processed batch ---

# The water removed is equal to the moisture content percentage of the fresh fruit weight
# Weight_fresh * moisture_content = water_removed_kg
weight_fresh_fruit_kg = water_removed_kg / moisture_content_fresh

# The final dry weight is the fresh weight minus the water removed
weight_dry_fruit_kg = weight_fresh_fruit_kg * (1 - moisture_content_fresh)

# --- Step 3: Calculate the total Sustainable Value Added (SVA) for the batch ---

# Total cost is the sum of costs for fruit, water, and energy
total_costs = cost_fruit + cost_water + cost_energy

# Total SVA for the batch is Revenue - Total Costs
total_sva_batch = total_revenue - total_costs

# --- Step 4: Calculate the SVA for the 250g product ---

# Convert the specific product's weight to kilograms
product_weight_kg = product_weight_g / 1000.0

# Calculate the SVA per kg of dry fruit
sva_per_kg = total_sva_batch / weight_dry_fruit_kg

# Calculate the final SVA for the 250g product
sva_for_product = sva_per_kg * product_weight_kg

# --- Step 5: Print the results and the final equation ---

print("Calculation of Sustainable Value Added (SVA) for a 250g Product\n")
print("Step 1: Calculate Total Batch Value Added")
print(f"  Total Revenue: ${total_revenue:.2f}")
print(f"  Total Costs (Fruit + Water + Energy): ${cost_fruit:.2f} + ${cost_water:.2f} + ${cost_energy:.2f} = ${total_costs:.2f}")
print(f"  Total SVA for Batch: ${total_revenue:.2f} - ${total_costs:.2f} = ${total_sva_batch:.2f}\n")

print("Step 2: Calculate Total Batch Weight")
print(f"  Water Removed: {water_removed_kg:.1f} kg")
print(f"  Fresh Fruit Moisture: {moisture_content_fresh*100:.0f}%")
print(f"  Total Dry Batch Weight = ({water_removed_kg:.1f} kg / {moisture_content_fresh}) * (1 - {moisture_content_fresh}) = {weight_dry_fruit_kg:.1f} kg\n")

print("Step 3: Calculate SVA for the 250g Product")
print("  SVA = (Total Batch SVA / Total Batch Weight) * Product Weight")
print(f"  SVA = (${total_sva_batch:.2f} / {weight_dry_fruit_kg:.1f} kg) * {product_weight_kg:.3f} kg")
print(f"  Final SVA = ${sva_for_product:.6f}")


# Restore original stdout
sys.stdout = old_stdout
# Print the captured output to the actual console
print(captured_output.getvalue())
# Final Answer format
print(f'<<<{sva_for_product:.6f}>>>')