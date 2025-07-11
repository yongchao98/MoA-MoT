import sys

# Given Parameters
travel_speed_mm_per_min = 30.0  # mm/min
arc_gap_mm = 6.0
# TIG welding arc efficiency (eta) is typically 0.6 to 0.8. We'll use 0.7 for calculation.
arc_efficiency = 0.7

# Answer Choices (Current in Amps, Voltage in Volts)
choices = {
    'A': {'I': 100.0, 'V': 7.5},
    'B': {'I': 17.5, 'V': 7.5},
    'C': {'I': 100.0, 'V': 15.0},
    'D': {'I': 150.0, 'V': 3.0},
    'E': {'I': 25.0, 'V': 17.5},
    'F': {'I': 80.0, 'V': 10.0}
}

# --- Step 1: Analyze Voltage based on Arc Gap ---
# A large 6 mm arc gap requires a relatively high voltage to sustain the arc.
# A common estimation is V = 10 + Arc_Gap_mm.
estimated_voltage = 10 + arc_gap_mm
print(f"Analysis Step 1: Evaluating Voltage")
print(f"For a {arc_gap_mm} mm arc gap, the estimated voltage is around {estimated_voltage:.1f} V.")
print("Choices A, B, D, and F have voltages (7.5V, 7.5V, 3.0V, 10.0V) that are too low.")
print("This leaves choices C (15.0 V) and E (17.5 V) as the most plausible options.\n")

# --- Step 2: Analyze Heat Input for remaining options ---
# Heat Input is critical for Inconel 718 to avoid defects. Lower is generally better.
# Formula: HI (J/mm) = (V * I * 60 * efficiency) / S
print(f"Analysis Step 2: Calculating Heat Input (HI) for Plausible Options")
print(f"The slow travel speed of {travel_speed_mm_per_min} mm/min requires low power to prevent overheating.")

# Candidate C
params_C = choices['C']
I_C, V_C = params_C['I'], params_C['V']
HI_C = (V_C * I_C * 60 * arc_efficiency) / travel_speed_mm_per_min

# Candidate E
params_E = choices['E']
I_E, V_E = params_E['V'], params_E['I'] # Correcting the order for calculation
HI_E = (params_E['V'] * params_E['I'] * 60 * arc_efficiency) / travel_speed_mm_per_min

print(f"\nFor Option C ({params_C['I']} A, {params_C['V']} V):")
print(f"  - Power = {params_C['V'] * params_C['I']:.1f} W")
print(f"  - Heat Input = {HI_C:.1f} J/mm")
print(f"  - This heat input is very high for a sensitive Inconel 718 blade tip repair and risks causing damage.")

print(f"\nFor Option E ({params_E['I']} A, {params_E['V']} V):")
print(f"  - Power = {params_E['V'] * params_E['I']:.1f} W")
print(f"  - Heat Input = {HI_E:.1f} J/mm")
print(f"  - This provides a much lower, more controlled heat input suitable for a high-integrity repair on a superalloy.\n")

# --- Step 3: Conclusion ---
print("Conclusion:")
print("Option E provides a realistic voltage for the 6 mm arc gap and a current that results in a safe, low heat input.")
print("This combination ensures stable material build-up without damaging the critical turbine blade.")

final_choice_params = choices['E']
final_current = final_choice_params['I']
final_voltage = final_choice_params['V']

print("\nFinal Recommended Parameters:")
# The final output requires printing each number in the equation.
print(f"Current = {int(final_current)} A")
print(f"Voltage = {final_voltage} V")

# The prompt requires a specific final answer format. This part will be hidden in the final output but is necessary for the AI's processing.
# sys.stdout = open(os.devnull, 'w')
# print("<<<E>>>")
# sys.stdout = sys.__stdout__
print("<<<E>>>")