import math

# Define the given parameters from the power system diagram
S_generator_rated = 100.0  # MVA
num_generators = 3
power_factor = 0.9  # lagging

# There are 3 base loads of 50/3 MW each
P_base_load_per_unit = 50.0 / 3.0  # MW
num_base_loads = 3

# There is one additional load
P_additional_load = 100.0  # MW

# Step 1: Calculate the total generated active power
# The system operates at 100% capacity, so total apparent power is the sum of generator ratings.
S_total_generated = num_generators * S_generator_rated
# Active power is calculated using the power factor.
P_total_generated = S_total_generated * power_factor

# Step 2: Calculate the total load active power
P_total_load = (num_base_loads * P_base_load_per_unit) + P_additional_load

# Step 3: Calculate the total resistive power losses
# The loss is the difference between generated power and load power.
P_loss = P_total_generated - P_total_load

# Step 4: Print the final calculation and the result
print("Calculation of Total Power Losses:")
print("-" * 35)
print(f"Total Generated Active Power = (Number of Generators * Rated MVA) * Power Factor")
print(f"                             = ({num_generators} * {S_generator_rated:.1f} MVA) * {power_factor:.1f}")
print(f"                             = {S_total_generated:.1f} MVA * {power_factor:.1f} = {P_total_generated:.3f} MW\n")

print(f"Total Load Active Power = (Number of Base Loads * Base Load Power) + Additional Load")
print(f"                        = ({num_base_loads} * {P_base_load_per_unit:.3f} MW) + {P_additional_load:.1f} MW")
print(f"                        = {num_base_loads * P_base_load_per_unit:.1f} MW + {P_additional_load:.1f} MW = {P_total_load:.3f} MW\n")

print(f"Total Power Loss = Total Generated Power - Total Load Power")
print(f"Final Equation: P_loss = {P_total_generated:.3f} MW - {P_total_load:.3f} MW")
print(f"Result: P_loss = {P_loss:.3f} MW")
print("-" * 35)