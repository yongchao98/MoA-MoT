import math

# Step 1: Define constants and given parameters
Io = 1e-9  # Reverse saturation current in Amperes
n = 1.5    # Diode ideality factor
T = 300    # Ambient temperature in Kelvin
V1 = 0.78  # Start voltage of the linear region in Volts
V2 = 0.98  # End voltage of the linear region in Volts
I2 = 0.445 # Current at V2 in Amperes
RL = 50.0  # Load resistance in Ohms
margin = 0.20 # 20% startup margin

# Physical constants
k = 1.380649e-23  # Boltzmann constant in J/K
q = 1.60217663e-19 # Elementary charge in Coulombs

# Step 2: Calculate the thermal voltage (Vt)
Vt = (k * T) / q

# Step 3: Calculate the current I1 at voltage V1 using the Shockley diode equation
# I1 = Io * (exp(V1 / (n * Vt)) - 1)
# The '-1' is negligible as the exponential term is very large
I1 = Io * math.exp(V1 / (n * Vt))

# Step 4: Calculate the dynamic resistance (rd) in the linear region
# rd = dV / dI = (V2 - V1) / (I2 - I1)
delta_V = V2 - V1
delta_I = I2 - I1
rd = delta_V / delta_I

# Step 5: For max power transfer from a negative resistance source, the load should be |rd|.
# Apply the 20% startup margin. For a negative resistance source, this means
# the transformed load should be less than |rd| to ensure startup.
target_impedance = abs(rd) * (1 - margin)

# Step 6: Calculate the impedance transformation ratio
# Ratio = Z_diode_side / Z_load_side
transformation_ratio = target_impedance / RL

# --- Output the results ---
print("--- Calculation Steps ---")
print(f"1. Thermal Voltage (Vt): {Vt:.4f} V")
print(f"2. Calculated Current at V1 (I1): {I1:.4f} A")
print(f"3. Dynamic Resistance (rd) = ({V2:.2f} V - {V1:.2f} V) / ({I2:.3f} A - {I1:.4f} A) = {rd:.4f} Ohms")
print(f"4. Target Impedance (with {margin*100}% margin) = |{rd:.4f}| * (1 - {margin:.2f}) = {target_impedance:.4f} Ohms")
print("\n--- Final Equation ---")
print(f"Transformation Ratio = Target Impedance / Load Resistance")
print(f"Transformation Ratio = {target_impedance:.4f} Ohms / {RL:.1f} Ohms")
print(f"\nFinal Answer: {transformation_ratio:.5f}")

# The final answer in the required format
final_answer = transformation_ratio
print(f"\n<<<{final_answer:.5f}>>>")