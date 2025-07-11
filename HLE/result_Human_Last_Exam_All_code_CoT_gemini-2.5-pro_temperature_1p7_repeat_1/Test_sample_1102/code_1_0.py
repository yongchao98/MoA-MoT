#
# Task: Determine the correct TIG welding parameters for a turbine blade repair.
#

# --- Step 1: Analyze the problem parameters ---
# Material: Inconel 718 (thin aerofoil section)
# Process: Manual TIG (GTAW)
# Arc Gap: 6 mm
# Travel Speed: 0.5 mm/s (slow)

# --- Step 2: Logical Deduction ---

# Deduction about Voltage based on Arc Gap:
# A large arc gap of 6 mm requires a relatively high voltage to maintain a stable electric arc.
# Voltages below ~15V are unlikely to sustain such a long arc.
# This points towards choices with higher voltages.
print("Analysis Step 1: The large 6 mm arc gap requires a relatively high voltage (likely > 15 V) to maintain a stable arc.")

# Deduction about Current based on material, geometry, and speed:
# Inconel 718 is heat-sensitive. The blade is thin, making it prone to burn-through and warping.
# The travel speed of 0.5 mm/s is very slow, which increases heat input per unit length.
# To prevent damage, a low heat input is critical, which means the welding current must be low.
print("Analysis Step 2: The thin Inconel 718 material and slow travel speed require a low current to prevent overheating and damage.")

# --- Step 3: Evaluate and Conclude ---
# We need a combination of LOW current and HIGH voltage.
# Reviewing the choices, 25 A and 17.5 V is the only option that fits both criteria.
chosen_current_amps = 25
chosen_voltage_volts = 17.5
chosen_option = "E"

print("\nConclusion: The best parameters combine low current with high voltage.")
print(f"The most suitable parameters are {chosen_current_amps} A and {chosen_voltage_volts} V.")

# --- Final Answer ---
# The final selected answer is E.
print("\nFinal selected equation values:")
print(f"Current = {chosen_current_amps} A")
print(f"Voltage = {chosen_voltage_volts} V")
print(f"<<<{chosen_option}>>>")
