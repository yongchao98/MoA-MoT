# --- Step 1: Define variables from the problem statement ---

# Economic and product-specific values
selling_price = 50.0  # in $
product_weight_kg = 0.250  # 250g converted to kg

# Actual resource usage for the 250g product
actual_water_L = 10.0
actual_energy_kWh = 15.0

# Benchmark (optimal) process details
# Assumption: Benchmark figures are for a standard 1kg of product
benchmark_output_kg = 1.0
benchmark_water_L = 100.0
benchmark_energy_kWh = 40.0
benchmark_water_cost = 10.0  # Cost for 100L of water
benchmark_energy_cost = 15.0  # Cost for 40kWh of energy

# --- Step 2: Calculate the price per unit for each resource ---
price_per_L_water = benchmark_water_cost / benchmark_water_L
price_per_kWh_energy = benchmark_energy_cost / benchmark_energy_kWh

# --- Step 3: Calculate the expected (benchmark) resource use for the 250g product ---
expected_water_L = (benchmark_water_L / benchmark_output_kg) * product_weight_kg
expected_energy_kWh = (benchmark_energy_kWh / benchmark_output_kg) * product_weight_kg

# --- Step 4: Calculate the total environmental cost impact ---
# This is the monetized value of the difference between actual and expected use.
# A negative value represents an environmental surplus (a good thing).
environmental_cost_impact = ((actual_water_L - expected_water_L) * price_per_L_water) + \
                             ((actual_energy_kWh - expected_energy_kWh) * price_per_kWh_energy)

# --- Step 5: Calculate the Sustainable Value Added (SVA) ---
# SVA = Economic Value Added - Environmental Cost Impact
sva = selling_price - environmental_cost_impact

# --- Final Output: Print the equation with all the numbers ---
print("The sustainable value added (SVA) is calculated as follows:")
print("SVA = Selling Price - Environmental Cost Impact")
print("Where Environmental Cost Impact = ((Actual Water - Expected Water) * Water Price) + ((Actual Energy - Expected Energy) * Energy Price)")
print("\nFinal Equation:")
print(f"SVA = {selling_price} - (({actual_water_L} - {expected_water_L}) * {price_per_L_water} + ({actual_energy_kWh} - {expected_energy_kWh}) * {price_per_kWh_energy})")
print(f"SVA = {selling_price} - ({environmental_cost_impact:.3f})")
print(f"SVA = {sva:.3f}")