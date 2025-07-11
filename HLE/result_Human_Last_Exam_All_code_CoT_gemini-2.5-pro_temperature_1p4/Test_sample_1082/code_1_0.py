# Step 1: Define the key financial values from the problem statement.
revenue = 50
cost_water = 10
cost_energy = 15

# Step 2: Calculate the total cost of the environmental resources.
total_resource_cost = cost_water + cost_energy

# Step 3: Calculate the Sustainable Value Added (SVA).
sva = revenue - total_resource_cost

# Step 4: Print the final equation, showing each number used in the calculation.
print("Sustainable Value Added (SVA) is calculated as Revenue minus the sum of resource costs.")
print("The final equation is:")
print(f"{revenue} - ({cost_water} + {cost_energy}) = {sva}")