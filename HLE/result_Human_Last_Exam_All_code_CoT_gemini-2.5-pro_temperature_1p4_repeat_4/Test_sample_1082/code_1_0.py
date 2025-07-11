# --- Define the given values ---

# Value Added (Economic Benefit)
selling_price = 50.0  # in dollars

# Corporate Eco-efficiency Benchmark (Optimal Use)
benchmark_water_volume = 100.0  # in liters
benchmark_water_cost = 10.0     # in dollars
benchmark_energy_volume = 40.0  # in kWh
benchmark_energy_cost = 15.0    # in dollars

# Actual Resource Consumption (Interpreted from the prompt)
actual_water_used = 10.0  # in liters
actual_energy_used = 15.0 # in kWh

# Note: The following variables are provided in the prompt but are not needed for this calculation
# as the resource usage is given in absolute terms, not as a rate per weight.
# net_weight_g = 250
# moisture_content = 0.20
# mango_ratio = 1
# coconut_ratio = 2

# --- Step 1: Calculate the unit cost of each resource from the benchmark ---
cost_per_liter_water = benchmark_water_cost / benchmark_water_volume
cost_per_kwh_energy = benchmark_energy_cost / benchmark_energy_volume

# --- Step 2: Calculate the total cost of the actual resources consumed ---
environmental_cost_water = actual_water_used * cost_per_liter_water
environmental_cost_energy = actual_energy_used * cost_per_kwh_energy
total_environmental_cost = environmental_cost_water + environmental_cost_energy

# --- Step 3: Calculate the Sustainable Value Added (SVA) ---
sva = selling_price - total_environmental_cost

# --- Step 4: Output the results showing the full equation ---
print("Sustainable Value Added (SVA) Calculation")
print("SVA = Value Added - Environmental Cost")
print("SVA = Selling Price - (Cost of Water Used + Cost of Energy Used)")
print(f"SVA = ${selling_price:.2f} - (${environmental_cost_water:.2f} + ${environmental_cost_energy:.3f})")
print(f"SVA = ${selling_price:.2f} - ${total_environmental_cost:.3f}")
print(f"Final Sustainable Value Added (SVA) = ${sva:.3f}")