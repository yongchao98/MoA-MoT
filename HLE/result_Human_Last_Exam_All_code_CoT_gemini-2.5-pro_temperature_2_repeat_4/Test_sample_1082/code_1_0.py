import math

# Step 1: Define initial parameters from the problem description
net_weight_g = 250.0
moisture_content = 0.20
mango_ratio = 1.0
coconut_ratio = 2.0
total_ratio = mango_ratio + coconut_ratio

selling_price = 50.0

# Resource consumption rates (interpreting "cost 10 and 10 and 15")
water_consumption_rate_L_per_kg = 10.0 # L per kg of total fresh fruit
mango_energy_rate_kWh_per_kg = 10.0      # kWh per kg of fresh mango
coconut_energy_rate_kWh_per_kg = 15.0    # kWh per kg of fresh coconut

# Benchmark resource usage for creating value
benchmark_water_L = 100.0
benchmark_energy_kWh = 40.0

# Step 2: Calculate the required weight of fresh fruits
net_weight_kg = net_weight_g / 1000.0
total_fresh_weight_kg = net_weight_kg / (1 - moisture_content)

fresh_mango_kg = total_fresh_weight_kg * (mango_ratio / total_ratio)
fresh_coconut_kg = total_fresh_weight_kg * (coconut_ratio / total_ratio)

# Step 3: Calculate the actual resources consumed for the product
actual_water_L = water_consumption_rate_L_per_kg * total_fresh_weight_kg
actual_energy_mango_kWh = mango_energy_rate_kWh_per_kg * fresh_mango_kg
actual_energy_coconut_kWh = coconut_energy_rate_kWh_per_kg * fresh_coconut_kg
total_actual_energy_kWh = actual_energy_mango_kWh + actual_energy_coconut_kWh

# Step 4: Calculate the benchmark value per unit of resource (shadow price)
benchmark_value_per_water = selling_price / benchmark_water_L
benchmark_value_per_energy = selling_price / benchmark_energy_kWh

# Step 5: Calculate the monetized environmental cost of the actual process
environmental_cost_water = actual_water_L * benchmark_value_per_water
environmental_cost_energy = total_actual_energy_kWh * benchmark_value_per_energy
total_environmental_cost = environmental_cost_water + environmental_cost_energy

# Step 6: Calculate the final Sustainable Value Added (SVA)
sva = selling_price - total_environmental_cost

# Print the final equation with each number as requested
print("The Sustainable Value Added (SVA) is calculated as follows:")
print("SVA = Selling Price - (Environmental Cost of Water + Environmental Cost of Energy)")
print(f"SVA = ${selling_price:.2f} - (${environmental_cost_water:.2f} + ${environmental_cost_energy:.2f})")
print(f"SVA = ${sva:.2f}")
print("\nFinal Calculation:")
# The final answer format requests the equation with the values substituted.
print(f"{sva:.2f} = {selling_price:.2f} - ({environmental_cost_water:.2f} + {environmental_cost_energy:.2f})")

<<<43.23>>>