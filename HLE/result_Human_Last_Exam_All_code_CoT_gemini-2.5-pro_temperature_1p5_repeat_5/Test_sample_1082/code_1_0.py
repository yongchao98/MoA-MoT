# Step 1: Define the key financial and cost variables from the problem.
revenue = 50.0
cost_water = 10.0
cost_energy = 15.0

# Step 2: Calculate the total environmental cost.
# This is the sum of the costs of all environmental resources used.
total_environmental_cost = cost_water + cost_energy

# Step 3: Calculate the Sustainable Value Added (SVA).
# SVA = Revenue - Total Environmental Cost
sva = revenue - total_environmental_cost

# Step 4: Print the breakdown of the calculation.
# The output clearly shows each number used in the final equation.
print("Calculating the Sustainable Value Added (SVA):")
print(f"Revenue from product = ${revenue}")
print(f"Cost of Water = ${cost_water}")
print(f"Cost of Energy = ${cost_energy}")
print(f"Total Environmental Cost = ${cost_water} + ${cost_energy} = ${total_environmental_cost}")
print(f"Final Equation: SVA = Revenue - Total Environmental Cost")
print(f"SVA = ${revenue} - ${total_environmental_cost}")
print(f"The Sustainable Value Added is ${sva}")
