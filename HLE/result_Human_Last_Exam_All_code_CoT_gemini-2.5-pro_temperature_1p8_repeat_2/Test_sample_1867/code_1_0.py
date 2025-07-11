import math

# Define the given parameters and constants
Io = 1e-9  # Reverse saturation current in Amperes
n = 1.5  # Diode ideality factor
T = 300  # Ambient temperature in Kelvin
V1 = 0.78  # Start voltage of linear region in Volts
V2 = 0.98  # End voltage of linear region in Volts
I2 = 0.445  # Current at V2 in Amperes
R_load_val = 50.0  # Load resistance in Ohms
margin = 0.20  # Startup margin of 20%

# Physical constants
k = 1.380649e-23  # Boltzmann constant in J/K
q = 1.60217663e-19  # Elementary charge in Coulombs

# --- Step 1: Calculate the initial state of the linear region ---
# Calculate thermal voltage (Vt)
Vt = (k * T) / q
print(f"Calculated Thermal Voltage (Vt): {Vt:.5f} V")

# Calculate the diode current I1 at voltage V1 using the diode equation
# This is the starting point of the linear behavior.
nVt = n * Vt
I1 = Io * (math.exp(V1 / nVt) - 1)
print(f"Calculated Diode Current (I1) at V1={V1}V: {I1:.5f} A")

# --- Step 2: Calculate the dynamic source resistance of the diode ---
delta_V = V2 - V1
delta_I = I2 - I1
r_d = delta_V / delta_I
print(f"Change in Voltage (ΔV): {delta_V:.2f} V")
print(f"Change in Current (ΔI): {delta_I:.5f} A")
print(f"Diode Dynamic Source Resistance (r_d = ΔV/ΔI): {r_d:.4f} Ohms")

# --- Step 3: Determine the target load impedance for optimal, stable power transfer ---
# For a negative resistance source, the load must be larger than the magnitude
# of the source resistance for stability. The margin is applied here.
r_d_mag = abs(r_d)
R_target = r_d_mag * (1 + margin)
print("\nFor stable operation with a negative resistance source, the target load must be greater than the source's magnitude.")
print(f"Applying a {margin*100}% startup margin...")
print("Target Impedance (R_target) = |Source Resistance| * (1 + Margin)")
print(f"R_target = {r_d_mag:.4f} * (1 + {margin:.2f}) = {R_target:.4f} Ohms")

# --- Step 4: Calculate the required impedance transformation ratio ---
# The ratio is from the load to the diode (i.e., diode side / load side)
ratio = R_target / R_load_val
print("\nThe impedance transformation ratio is the ratio of the target impedance to the load impedance.")
print("Transformation Ratio = Target Impedance / Load Impedance")
print(f"Transformation Ratio = {R_target:.4f} / {R_load_val:.1f} = {ratio:.5f}")

# Final Answer
print(f"\nThe final required impedance transformation ratio is {ratio:.5f}.")
print(f'<<<{ratio:.5f}>>>')