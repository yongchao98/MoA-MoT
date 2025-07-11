import math

# --- Given Parameters ---
Io = 1e-9  # Reverse saturation current in Amperes
n = 1.5    # Diode ideality factor
T = 300    # Ambient temperature in Kelvin
V1 = 0.78  # Start voltage of linear region in Volts
V2 = 0.98  # End voltage of linear region in Volts
I2 = 0.445 # Current at V2 in Amperes
RL = 50.0  # Load resistance in Ohms
margin = 0.20 # 20% startup margin

# --- Physical Constants ---
k = 1.380649e-23  # Boltzmann constant in J/K
q = 1.60217663e-19  # Elementary charge in Coulombs

# --- Step 1: Calculate Thermal Voltage (Vt) ---
Vt = (k * T) / q

# --- Step 2: Calculate Current I1 at V1 ---
# Assuming normal diode behavior up to V1, we use the Shockley equation:
# I = Io * (exp(V / (n * Vt)) - 1)
exponent = V1 / (n * Vt)
I1 = Io * (math.exp(exponent) - 1)

# --- Step 3: Calculate Dynamic Resistance (rd) ---
# In the linear region from (V1, I1) to (V2, I2), the resistance is the slope dV/dI.
delta_V = V2 - V1
delta_I = I2 - I1
rd = delta_V / delta_I
abs_rd = abs(rd)

# --- Step 4: Calculate Target Input Impedance with Margin ---
# For oscillator startup, the load R_in must be > |rd|. We apply the margin.
R_in_target = abs_rd * (1 + margin)

# --- Step 5: Calculate Impedance Transformation Ratio ---
transformation_ratio = R_in_target / RL

# --- Output the results ---
print("This problem involves a diode with a negative dynamic resistance, used as a signal source.")
print("The goal is to find the impedance ratio for optimal power transfer with a startup margin.")

print("\n--- Calculation Breakdown ---")
print(f"1. Calculated Diode Current at V1={V1}V (I1): {I1:.5f} A")
print(f"2. Calculated Dynamic Resistance (rd): {rd:.5f} Ohms")
print(f"3. Calculated Target Input Impedance with 20% Margin (R_in): {R_in_target:.5f} Ohms")
print(f"4. Given Load Resistance (RL): {RL:.1f} Ohms")

print("\n--- Final Equation ---")
print("Impedance Transformation Ratio = (Target Input Impedance) / (Load Resistance)")
print(f"Impedance Transformation Ratio = {R_in_target:.5f} / {RL:.1f}")

print(f"\nFinal Answer: {transformation_ratio:.5f}")
print(f"<<<{transformation_ratio:.5f}>>>")