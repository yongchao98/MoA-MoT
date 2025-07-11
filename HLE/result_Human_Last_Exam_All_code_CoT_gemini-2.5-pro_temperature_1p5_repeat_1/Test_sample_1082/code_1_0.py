# Step 1: Define the given parameters
net_weight_dry_g = 250.0  # Net weight of the final dry product in grams
selling_price = 50.0  # Selling price of the product in dollars
moisture_content_fresh = 0.20  # Moisture content of the fresh fruits

# Benchmark values are assumed to be per 1 kg of FRESH fruit input
benchmark_water_use_L_per_kg_fresh = 100.0
benchmark_energy_use_kWh_per_kg_fresh = 40.0

# Interpreting "cost 10 and 10 and 15 water and energy respectively" as:
# $10 for the benchmark water usage (100L) and $15 for the benchmark energy usage (40 kWh).
benchmark_total_water_cost = 10.0
benchmark_total_energy_cost = 15.0

# --- Calculations ---

# Step 2: Calculate the mass of fresh fruit needed
# Convert dry weight from grams to kilograms
net_weight_dry_kg = net_weight_dry_g / 1000.0
# The solid content of the fresh fruit is (1 - moisture content)
solid_content_fraction = 1 - moisture_content_fresh
# Calculate the initial mass of fresh fruit required
fresh_fruit_mass_kg = net_weight_dry_kg / solid_content_fraction

# Step 3: Calculate actual resource consumption for the product
actual_water_use_L = fresh_fruit_mass_kg * benchmark_water_use_L_per_kg_fresh
actual_energy_use_kWh = fresh_fruit_mass_kg * benchmark_energy_use_kWh_per_kg_fresh

# Step 4: Calculate the total Environmental Cost (EC)
# First, find the cost per unit of each resource
cost_per_liter_water = benchmark_total_water_cost / benchmark_water_use_L_per_kg_fresh
cost_per_kWh_energy = benchmark_total_energy_cost / benchmark_energy_use_kWh_per_kg_fresh

# Then, calculate the cost for the actual resources used
environmental_cost_water = actual_water_use_L * cost_per_liter_water
environmental_cost_energy = actual_energy_use_kWh * cost_per_kWh_energy
total_environmental_cost = environmental_cost_water + environmental_cost_energy

# Step 5: Calculate the Sustainable Value Added (SVA)
# Value Added (VA) is the selling price
value_added = selling_price
sva = value_added - total_environmental_cost

# --- Final Output ---
print("Calculating Sustainable Value Added (SVA):")
print(f"Value Added (Selling Price) = ${value_added:.2f}")
print(f"Environmental Cost = ${total_environmental_cost:.2f}")
print("-" * 30)
# Final equation as requested, showing the main numbers
print(f"SVA = ${value_added:.2f} - ${total_environmental_cost:.2f} = ${sva:.2f}")
