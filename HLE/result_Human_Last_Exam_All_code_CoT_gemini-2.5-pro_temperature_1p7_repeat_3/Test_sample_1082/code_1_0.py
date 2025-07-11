import sys

# --- Step 1: Define initial values from the problem description ---

# Value Added is the revenue from the specific product
product_revenue = 50.0  # in dollars
product_net_weight_g = 250.0  # in grams

# Environmental costs for the entire benchmark batch process
total_water_cost_for_batch = 10.0  # in dollars
total_energy_cost_for_batch = 15.0  # in dollars

# Details of the dehydration process
water_removed_liters = 100.0  # in liters
fresh_fruit_moisture_percentage = 0.20  # 20% moisture content

# --- Step 2: Calculate the total mass of the dried fruit from the benchmark process ---

# Assuming density of water is 1 kg/L, so 100L of water is 100kg
water_removed_kg = water_removed_liters

# If the removed water is 20% of the total fresh mass
total_fresh_mass_kg = water_removed_kg / fresh_fruit_moisture_percentage

# The final dry mass is the initial mass minus the water removed
total_dry_mass_kg = total_fresh_mass_kg - water_removed_kg

# --- Step 3: Calculate the total environmental cost for the batch ---
total_environmental_cost = total_water_cost_for_batch + total_energy_cost_for_batch

# --- Step 4: Calculate the environmental cost for the specific 250g product ---

# Convert product weight to kg to match batch units
product_net_weight_kg = product_net_weight_g / 1000.0

# Calculate the environmental cost proportionally
product_environmental_cost = (total_environmental_cost / total_dry_mass_kg) * product_net_weight_kg

# --- Step 5: Calculate the final Sustainable Value Added (SVA) ---

# SVA = Value Added - Environmental Cost
# Here, we use product_revenue as the Value Added
sustainable_value_added = product_revenue - product_environmental_cost

# --- Final Output ---
# Print the final equation with all the calculated numbers
print("Sustainable Value Added (SVA) Calculation:")
print(f"Value Added = ${product_revenue:.2f}")
print(f"Environmental Cost = ${product_environmental_cost:.6f}")
print(f"SVA = Value Added - Environmental Cost")
sys.stdout.write(f"SVA = ${product_revenue:.2f} - ${product_environmental_cost:.6f} = ${sustainable_value_added:.6f}\n")
sys.stdout.flush()
# The ratios of mangoes to coconuts, and specific energy/water usage numbers (40 kWh, 100L)
# were not needed for the final calculation because the total environmental cost ($10 + $15)
# for the entire batch was provided directly.