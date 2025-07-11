import math

# --- Given Parameters ---

# Generator parameters
s_gen_per_unit = 100.0  # MVA per generator
num_generators = 3
power_factor = 0.9      # lagging

# Load parameters
p_load_base_total = 50.0  # MW (3 loads of 50/3 MW each)
p_load_additional = 100.0 # MW

# --- Calculation ---

# Step 1: Calculate total generated active power
total_s_gen = num_generators * s_gen_per_unit
total_p_gen = total_s_gen * power_factor

# Step 2: Calculate total consumed active power by loads
total_p_load = p_load_base_total + p_load_additional

# Step 3: Calculate the total power loss
total_p_loss = total_p_gen - total_p_load

# --- Output Results ---

print("This script calculates the total power loss using the power balance method.")

print("\n--- Calculation Steps ---")

print("\n1. Calculate Total Generated Active Power (P_generated):")
print(f"The system has {num_generators} generators, each operating at {s_gen_per_unit} MVA with a power factor of {power_factor}.")
print(f"Equation: P_generated = (Number of Generators * MVA per Generator) * Power Factor")
print(f"Calculation: ({num_generators} * {s_gen_per_unit:.1f}) MVA * {power_factor} = {total_p_gen:.1f} MW")

print("\n2. Calculate Total Consumed Active Power (P_load):")
print(f"The total base load is {p_load_base_total:.1f} MW and the additional load is {p_load_additional:.1f} MW.")
print(f"Equation: P_load = Total Base Load + Additional Load")
print(f"Calculation: {p_load_base_total:.1f} MW + {p_load_additional:.1f} MW = {total_p_load:.1f} MW")

print("\n3. Calculate Total Resistive Power Loss (P_loss):")
print(f"Equation: P_loss = P_generated - P_load")
print(f"Final Calculation: {total_p_gen:.1f} MW - {total_p_load:.1f} MW = {total_p_loss:.3f} MW")

print("\n-------------------------------------------------")
print(f"The total power loss in the system is {total_p_loss:.3f} MW.")
print("-------------------------------------------------")
