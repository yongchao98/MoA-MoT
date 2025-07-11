import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Define Parameters from the Problem Description ---

# Product characteristics (context for the calculation)
net_weight_g = 250
mango_ratio = 1
coconut_ratio = 2
moisture_content_fresh = 0.20

# Financial and Environmental Impact Data
selling_price = 50
water_usage_liters = 100
energy_usage_kwh = 40
water_cost_factor = 10
energy_cost_factor = 15

# --- Step-by-Step Calculation ---

# Step 1: Calculate the environmental cost attributed to water usage.
water_environmental_cost = water_usage_liters * water_cost_factor

# Step 2: Calculate the environmental cost attributed to energy usage.
energy_environmental_cost = energy_usage_kwh * energy_cost_factor

# Step 3: Sum the individual environmental costs to get the total.
total_environmental_cost = water_environmental_cost + energy_environmental_cost

# Step 4: Calculate the final Sustainable Value Added (SVA).
sva = selling_price - total_environmental_cost

# --- Output the Results ---

# Print the calculation steps clearly as requested.
print("Calculating the Sustainable Value Added (SVA):")
print(f"SVA = Selling Price - Total Environmental Cost")
print(f"SVA = {selling_price} - (Water Usage * Water Cost Factor + Energy Usage * Energy Cost Factor)")
print(f"SVA = {selling_price} - ({water_usage_liters} * {water_cost_factor} + {energy_usage_kwh} * {energy_cost_factor})")
print(f"SVA = {selling_price} - ({water_environmental_cost} + {energy_environmental_cost})")
print(f"SVA = {selling_price} - {total_environmental_cost}")
print(f"The final Sustainable Value Added (SVA) is: {sva}")

# --- Final Answer Formatting ---
# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the captured output to the actual console
print(output)
# Print the final answer in the specified format
print(f"<<<{sva}>>>")