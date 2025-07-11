import math

# Given welding parameters
material = "Inconel 718"
process = "Manual TIG (GTAW) build-up"
travel_speed_mms = 0.5  # mm per second
arc_gap_mm = 6

# --- Analysis ---
# 1. Voltage Estimation: For a TIG arc in Argon, V is approx. 10-12V + (1V * arc_length_mm).
# For a 6mm arc gap, Voltage should be around 10 + 6 = 16V.
# From the choices, 15V and 17.5V are plausible. 3V, 7.5V and 10V are too low for a 6mm arc.

# 2. Current Estimation: For building up an Inconel 718 blade tip, a very low current (<50A)
# would be insufficient to create a stable weld pool. A very high current (>120A) would
# risk overheating and burn-through on a thin edge. Currents of 80A or 100A are plausible.

# 3. Combining the analysis:
#   - Choice A (100 A, 7.5 V): Voltage too low.
#   - Choice B (17.5 A, 7.5 V): Voltage and Current too low.
#   - Choice C (100 A, 15 V): Plausible voltage and current.
#   - Choice D (150 A, 3 V): Voltage too low.
#   - Choice E (25 A, 17.5 V): Plausible voltage, but current is too low for build-up.
#   - Choice F (80 A, 10 V): Voltage too low for a 6mm arc gap.

# Conclusion: Choice C is the most technically sound option.
selected_current_A = 100  # Amps
selected_voltage_V = 15   # Volts

# Calculate Heat Input (Joules per mm)
# Heat Input = (Voltage * Current * Efficiency) / Travel Speed
# We will assume an efficiency of 1 (or 100%) for the purpose of the equation,
# though real-world TIG efficiency is typically 0.6-0.8.
efficiency = 1.0
heat_input_J_per_mm = (selected_voltage_V * selected_current_A * efficiency) / travel_speed_mms

# --- Output ---
print("Welding Procedure Specification (WPS) Recommendation")
print("-" * 50)
print(f"Component: Rolls-Royce Trent 1000 Compressor Blade Tip")
print(f"Material: {material}")
print(f"Process: {process}")
print("-" * 50)
print("Based on the specified 6 mm arc gap, a voltage of approximately 15-17V is required to maintain a stable arc.")
print("For a root pass build-up on Inconel 718, a current of 100A provides sufficient energy for proper fusion without being excessive.")
print("\nRecommended Parameters for the Root Pass:")
print(f"Current: {selected_current_A} A")
print(f"Voltage: {selected_voltage_V} V")
print("\nThis corresponds to a calculated heat input based on the equation:")
print(f"Heat Input (J/mm) = ({selected_voltage_V} V * {selected_current_A} A) / {travel_speed_mms} mm/s = {heat_input_J_per_mm:.0f} J/mm")
print("-" * 50)
