import math

# --- Given Parameters ---
Io = 1e-9  # Reverse saturation current in Amperes
n = 1.5    # Diode ideality factor
T = 300    # Ambient temperature in Kelvin
V1 = 0.78  # Start voltage of the linear region in Volts
V2 = 0.98  # End voltage of the linear region in Volts
I2 = 0.445 # Current at V2 in Amperes
R_L = 50   # Load resistance in Ohms
margin = 0.20 # Startup margin

# --- Physical Constants ---
k = 1.380649e-23  # Boltzmann's constant in J/K
q = 1.60217663e-19  # Elementary charge in Coulombs

# Step 1: Calculate the thermal voltage (V_T)
V_T = (k * T) / q

# Step 2: Calculate the current I1 at V1 using the diode equation
# The exponential term is much larger than 1, so we can approximate I = Io * exp(...)
I1 = Io * (math.exp(V1 / (n * V_T)) - 1)

# Step 3: Calculate the dynamic resistance (rd) in the linear region
# rd = dV / dI = (V2 - V1) / (I2 - I1)
delta_V = V2 - V1
delta_I = I2 - I1
rd = delta_V / delta_I
abs_rd = abs(rd)

# Step 4: Apply the 20% startup margin to the dynamic resistance
# For startup, R_load_seen < |rd|. The margin makes the target impedance smaller.
R_target = abs_rd * (1 - margin)

# Step 5: Calculate the final impedance transformation ratio
# Ratio = Z_primary / Z_secondary = R_target / R_L
transformation_ratio = R_target / R_L

# --- Output the results ---
print(f"Step 1: The thermal voltage V_T is {V_T:.4f} V.")
print(f"Step 2: The current I1 at V1={V1}V is calculated to be {I1:.4f} A.")
print(f"Step 3: The dynamic resistance rd = ({V2}V - {V1}V) / ({I2}A - {I1:.4f}A) = {rd:.4f} Ohms.")
print(f"Step 4: The target impedance for matching, including a {margin*100}% startup margin, is |{rd:.4f}| * (1 - {margin}) = {R_target:.4f} Ohms.")
print("\nFinal Calculation:")
print(f"The required impedance transformation ratio is the target impedance divided by the load resistance.")
print(f"Ratio = {R_target:.4f} Ohms / {R_L} Ohms")
print(f"Result: {transformation_ratio:.5f}")
