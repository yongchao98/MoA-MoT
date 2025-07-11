# Variables from the problem
net_weight_g = 250.0
moisture_content = 0.20
revenue = 50.0

# Benchmark resource usage (ASSUMPTION: for 1kg or 1000g of fresh fruit)
benchmark_fresh_fruit_g = 1000.0
benchmark_water_L = 100.0
benchmark_energy_kWh = 40.0

# Environmental burden 'cost' factors
cost_factor_water = 10.0
cost_factor_energy = 15.0

# Step 1: Calculate the weight of fresh fruit needed for the 250g product.
# Dry weight is (1 - moisture_content) of the fresh weight.
# Fresh weight = Dry weight / (1 - moisture_content)
fresh_fruit_needed_g = net_weight_g / (1 - moisture_content)

# Step 2: Calculate resource consumption rates per gram of fresh fruit from the benchmark.
water_rate_L_per_g = benchmark_water_L / benchmark_fresh_fruit_g
energy_rate_kWh_per_g = benchmark_energy_kWh / benchmark_fresh_fruit_g

# Step 3: Calculate the total resources used for the product.
water_used_L = fresh_fruit_needed_g * water_rate_L_per_g
energy_used_kWh = fresh_fruit_needed_g * energy_rate_kWh_per_g

# Step 4: Calculate the environmental burden for each resource and the total burden.
water_burden = water_used_L * cost_factor_water
energy_burden = energy_used_kWh * cost_factor_energy
total_environmental_burden = water_burden + energy_burden

# Step 5: Calculate the final Sustainable Value Added (SVA).
sva = revenue - total_environmental_burden

# --- Output the results ---
print(f"1. Revenue from Product: ${revenue:.0f}")

print("\n2. Environmental Burden Calculation:")
print(f"  - Initial fresh fruit weight needed: {fresh_fruit_needed_g:.1f}g")
print(f"  - Total water used for product: {water_used_L:.2f}L")
print(f"  - Total energy used for product: {energy_used_kWh:.2f}kWh")
print(f"  - Water Burden = {water_used_L:.2f} * {cost_factor_water:.0f} = {water_burden:.1f}")
print(f"  - Energy Burden = {energy_used_kWh:.2f} * {cost_factor_energy:.0f} = {energy_burden:.1f}")
print(f"  - Total Environmental Burden = {water_burden:.1f} + {energy_burden:.1f} = {total_environmental_burden:.0f}")

print("\n3. Final Sustainable Value Added (SVA) Equation:")
print(f"SVA = Revenue - Total Environmental Burden")
print(f"{revenue:.0f} - {total_environmental_burden:.0f} = {sva:.0f}")