import math

# --- Given Parameters ---
Io = 1e-9  # Reverse saturation current in Amperes
n = 1.5    # Diode ideality factor
T = 300    # Ambient temperature in Kelvin
V1 = 0.78  # Start voltage of the linear region in Volts
V2 = 0.98  # End voltage of the linear region in Volts
I2 = 0.445 # Current at V2 in Amperes
R_load = 50 # Actual load resistance in Ohms
margin = 0.20 # Startup margin

# --- Physical Constants ---
k = 1.380649e-23 # Boltzmann's constant in J/K
q = 1.602176634e-19 # Elementary charge in Coulombs

# Step 1: Calculate thermal voltage (VT) and the current I1 at V1
Vt = (k * T) / q
# The problem implies the linear region starts from the normal diode behavior at V1.
# I = Io * (exp(V / (n * Vt)) - 1). The '-1' is negligible at this forward bias.
I1 = Io * (math.exp(V1 / (n * Vt)))

# Step 2: Calculate the dynamic resistance of the diode (the source)
# Dynamic resistance Rs = dV / dI.
# Since the region is defined as linear, we use delta V / delta I.
dV = V2 - V1
dI = I2 - I1
# For max power transfer, we match the load to the magnitude of the source's dynamic resistance.
Rs = abs(dV / dI)

# Step 3: Apply the 20% startup margin to find the target load impedance
R_target = Rs * (1 + margin)

# Step 4: Calculate the required impedance transformation ratio
# Ratio = Z_actual_load / Z_target_load
transformation_ratio = R_load / R_target

# --- Final Output ---
print("This script calculates the impedance transformation ratio for a diode signal source.")
print("\nStep 1: Calculate the starting current (I1) of the linear region.")
print(f"I1 = Io * exp(V1 / (n * Vt)) = {Io:.1e} * exp({V1} / ({n} * {Vt:.4f})) = {I1:.4f} A")

print("\nStep 2: Calculate the diode's dynamic source resistance (Rs).")
print(f"Rs = |(V2 - V1) / (I2 - I1)| = |({V2} - {V1}) / ({I2} - {I1:.4f})| = {Rs:.4f} Ohms")

print("\nStep 3: Apply the startup margin to determine the target impedance.")
print(f"R_target = Rs * (1 + margin) = {Rs:.4f} * (1 + {margin}) = {R_target:.4f} Ohms")

print("\nStep 4: Calculate the final transformation ratio.")
# The final equation showing all the numbers
print(f"Transformation Ratio = R_load / R_target")
print(f"Transformation Ratio = {R_load} / {R_target:.4f}")
print(f"Result: {transformation_ratio:.4f}")

# Final Answer in requested format
print(f"\n<<<final_answer>>>\nThe final equation is:\nTransformation Ratio = {R_load} / (|({V2} - {V1}) / ({I2} - {I1:.4f})| * (1 + {margin})) = {transformation_ratio:.4f}\n<<<{transformation_ratio:.2f}>>>")