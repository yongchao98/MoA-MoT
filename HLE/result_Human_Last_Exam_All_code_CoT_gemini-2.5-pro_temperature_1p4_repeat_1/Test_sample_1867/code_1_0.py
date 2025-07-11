import math

# Description:
# This script calculates the required impedance transformation ratio for a diode
# acting as a signal source, considering a startup margin for optimum power transfer.

# --- Given Parameters ---
Io = 1e-9       # A, reverse saturation current
n = 1.5         # diode ideality factor
T = 300.0       # K, ambient temperature
R_L = 50.0      # ohm, load resistance
V1 = 0.78       # V, start voltage of the linear region
V2 = 0.98       # V, end voltage of the linear region
I2 = 0.445      # A, current at V2
margin = 0.20   # 20% startup margin

# --- Physical Constants ---
k = 1.380649e-23  # J/K, Boltzmann constant
q = 1.60217663e-19 # C, elementary charge

# --- Step 1: Calculate the source's dynamic resistance (Rs) ---

# First, calculate the thermal voltage (Vt)
Vt = (k * T) / q

# Second, calculate the current I1 at voltage V1 using the diode equation
# I1 = Io * (exp(V1 / (n * Vt)) - 1). The '-1' is negligible but included for accuracy.
I1 = Io * (math.exp(V1 / (n * Vt)) - 1)

# Third, calculate the dynamic resistance (Rs) as the slope dV/dI
Rs = (V2 - V1) / (I2 - I1)
Rs_abs = abs(Rs)

print("Step 1: Calculate the source dynamic resistance (Rs).")
print(f"The calculation for the source resistance Rs = (V2 - V1) / (I2 - I1) requires I1.")
print(f"The current I1 at V1 = {V1} V is found using the diode equation:")
print(f"I1 = {Io:.1e} * (exp({V1} / ({n} * {Vt:.5f})) - 1) = {I1:.5f} A")
print(f"Now, Rs can be calculated:")
print(f"Rs = ({V2} - {V1}) / ({I2} - {I1:.5f}) = {Rs:.5f} Ω\n")


# --- Step 2: Determine the required transformed load resistance (RL_transformed) ---

# To ensure startup, |Rs| must be larger than the load it sees.
# With a 20% margin, we set: |Rs| = (1 + margin) * RL_transformed
RL_transformed = Rs_abs / (1 + margin)

print("Step 2: Determine the required transformed load resistance (RL_t) with a 20% margin.")
print("The condition for startup is |Rs| >= (1 + margin) * RL_t.")
print("We design based on this condition to ensure oscillation starts.")
print(f"RL_t = |Rs| / (1 + margin) = {Rs_abs:.5f} / (1 + {margin}) = {RL_transformed:.5f} Ω\n")


# --- Step 3: Calculate the final impedance transformation ratio ---

# The ratio is from the load to the diode: Z_diode / Z_load
transformation_ratio = RL_transformed / R_L

print("Step 3: Calculate the final impedance transformation ratio from load to diode.")
print("The ratio is the required impedance at the diode divided by the load impedance.")
print(f"Transformation Ratio = RL_t / R_L = {RL_transformed:.5f} Ω / {R_L} Ω = {transformation_ratio:.5f}")

<<<0.03384>>>