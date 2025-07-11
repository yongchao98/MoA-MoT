import math

# --- Problem Parameters ---

# Generation
NUM_GENERATORS = 3
S_PER_GENERATOR_MVA = 100
POWER_FACTOR = 0.9

# Loads
NUM_BASE_LOADS = 3
P_PER_BASE_LOAD_MW = 50 / 3
P_ADDITIONAL_LOAD_MW = 100

# --- Step-by-step Calculation ---

# Step 1: Calculate total active power load (P_Load)
p_total_base_load_mw = NUM_BASE_LOADS * P_PER_BASE_LOAD_MW
p_total_load_mw = p_total_base_load_mw + P_ADDITIONAL_LOAD_MW

# Step 2: Calculate total active power generated (P_Generated)
# "operates at 100% capacity" implies generators output their rated MVA.
s_total_generated_mva = NUM_GENERATORS * S_PER_GENERATOR_MVA
p_total_generated_mw = s_total_generated_mva * POWER_FACTOR

# Step 3: Calculate total resistive power loss (P_Loss) using the power balance equation
p_loss_mw = p_total_generated_mw - p_total_load_mw

# --- Output the results and the final equation ---
print("The total resistive power loss is calculated using the power balance equation: P_Loss = P_Generated - P_Load\n")

print(f"1. Calculate Total Generated Active Power (P_Generated):")
print(f"   Total Apparent Power = {NUM_GENERATORS} generators * {S_PER_GENERATOR_MVA} MVA = {s_total_generated_mva:.0f} MVA")
print(f"   P_Generated = Total Apparent Power * Power Factor")
print(f"   P_Generated = {s_total_generated_mva:.0f} MVA * {POWER_FACTOR} = {p_total_generated_mw:.0f} MW\n")

print(f"2. Calculate Total Active Power Load (P_Load):")
print(f"   Total Base Load = {NUM_BASE_LOADS} loads * {P_PER_BASE_LOAD_MW:.3f} MW = {p_total_base_load_mw:.0f} MW")
print(f"   P_Load = Total Base Load + Additional Load")
print(f"   P_Load = {p_total_base_load_mw:.0f} MW + {P_ADDITIONAL_LOAD_MW} MW = {p_total_load_mw:.0f} MW\n")

print("3. Calculate Total Power Loss (P_Loss):")
print("   P_Loss = P_Generated - P_Load")
print(f"   The final equation with the numbers is:")
print(f"   {p_total_generated_mw:.0f} MW - {p_total_load_mw:.0f} MW = {p_loss_mw:.3f} MW")

print(f"\nResult: The total power losses in the system are {p_loss_mw:.3f} MW.")
<<<120.000>>>