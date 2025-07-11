import sys
# Redirect stdout to a variable to prevent printing before the final format
original_stdout = sys.stdout
from io import StringIO
captured_output = StringIO()
sys.stdout = captured_output


# Plan:
# 1. Identify the economic value created from the sale of the product.
# 2. Identify the environmental costs (eco-costs) associated with the production process.
# 3. Calculate the Sustainable Value Added (SVA) by subtracting the total eco-cost from the economic value.
# 4. Print the final calculation, showing all the numbers involved in the equation.

# Step 1: Define the economic value.
# The final product is sold for $50.
value_created = 50.0

# Step 2: Define the eco-costs.
# The problem states the process costs $10 for water and $15 for energy.
eco_cost_water = 10.0
eco_cost_energy = 15.0

# Calculate the total eco-cost.
total_eco_cost = eco_cost_water + eco_cost_energy

# Step 3: Calculate the Sustainable Value Added (SVA).
sustainable_value_added = value_created - total_eco_cost

# Step 4: Print the final output showing the equation with all numbers.
# The format should clearly display how the SVA is derived.
print("Sustainable Value Added (SVA) Calculation:")
print(f"SVA = (Economic Value Created) - (Total Eco-Cost)")
print(f"SVA = ${value_created} - (${eco_cost_water} + ${eco_cost_energy})")
print(f"SVA = ${value_created} - ${total_eco_cost}")
print(f"The final Sustainable Value Added is ${sustainable_value_added}")

# Get the content from captured_output and restore stdout
output = captured_output.getvalue()
sys.stdout = original_stdout
print(output)
print(f'<<<${sustainable_value_added}>>>')
