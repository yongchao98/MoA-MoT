import sys

# --- Define initial values from the problem ---

# Product and Revenue
product_net_weight_g = 250
selling_price = 50.0

# Benchmark resource usage (for the entire process)
benchmark_water_L = 100.0
benchmark_energy_kWh = 40.0

# Resource Costs
cost_per_L_water = 10.0
cost_per_kWh_energy = 15.0

# --- State the assumption ---
# We assume the benchmark resource usage produces 1kg (1000g) of the final product.
benchmark_output_g = 1000.0

# --- Calculations ---

# 1. Calculate resource usage per gram of product
water_L_per_g = benchmark_water_L / benchmark_output_g
energy_kWh_per_g = benchmark_energy_kWh / benchmark_output_g

# 2. Calculate resource usage for the specific 250g product
water_for_product_L = water_L_per_g * product_net_weight_g
energy_for_product_kWh = energy_kWh_per_g * product_net_weight_g

# 3. Calculate the environmental cost for the 250g product
water_cost = water_for_product_L * cost_per_L_water
energy_cost = energy_for_product_kWh * cost_per_kWh_energy
total_environmental_cost = water_cost + energy_cost

# 4. Calculate the final Sustainable Value Added (SVA)
# SVA = Revenue - Total Environmental Cost
sva = selling_price - total_environmental_cost

# --- Output the results ---

print("The Sustainable Value Added (SVA) is calculated as Revenue - Environmental Costs.")
print(f"Revenue = ${selling_price:.2f}")
print(f"Environmental Costs = Water Cost + Energy Cost")
print(f"  - Water Cost = {water_for_product_L:.2f} L * ${cost_per_L_water:.2f}/L = ${water_cost:.2f}")
print(f"  - Energy Cost = {energy_for_product_kWh:.2f} kWh * ${cost_per_kWh_energy:.2f}/kWh = ${energy_cost:.2f}")

print("\nFinal Equation:")
# The user requested to output each number in the final equation.
# We redirect stdout to capture the string representation for the final format.
original_stdout = sys.stdout
from io import StringIO
captured_output = StringIO()
sys.stdout = captured_output

# This print statement is captured
print(f"SVA = {selling_price:.0f} - ({water_cost:.0f} + {energy_cost:.0f}) = {sva:.0f}")

# Restore stdout and print the captured string
sys.stdout = original_stdout
final_equation_str = captured_output.getvalue().strip()
print(final_equation_str)

# Also printing the final value for easy access
# print(f"\nFinal SVA = ${sva:.2f}")
