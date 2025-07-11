import math

# --- Step 1: Define constants and problem parameters ---
V1 = 0.78       # Start voltage in Volts
V2 = 0.98       # End voltage in Volts
I2 = 0.445      # Current at V2 in Amps
Io = 1e-9       # Reverse saturation current in Amps
n = 1.5         # Diode ideality factor
T = 300         # Ambient temperature in Kelvin
RL_load = 50    # Load resistance in Ohms
margin = 0.20   # Startup margin of 20%

# Physical constants
k = 1.380649e-23  # Boltzmann constant in J/K
q = 1.6021766e-19  # Elementary charge in Coulombs

# --- Step 2: Calculate the initial current I1 ---
# Calculate thermal voltage (Vt)
Vt = (k * T) / q
# Calculate current I1 at voltage V1 using the Shockley diode equation
I1 = Io * (math.exp(V1 / (n * Vt)) - 1)

print(f"Calculated thermal voltage (Vt): {Vt:.4f} V")
print(f"Calculated initial current (I1) at {V1} V: {I1:.4f} A")

# --- Step 3: Calculate the diode's dynamic resistance (Rs) ---
# Dynamic resistance is dV/dI in the linear region
delta_V = V2 - V1
delta_I = I2 - I1
Rs = delta_V / delta_I
Rs_mag = abs(Rs)

print(f"\nChange in Voltage (ΔV): {delta_V:.2f} V")
print(f"Change in Current (ΔI): {delta_I:.4f} A")
print(f"Diode dynamic resistance (Rs = ΔV/ΔI): {Rs:.4f} Ohms")
print(f"Magnitude of the diode's resistance |Rs|: {Rs_mag:.4f} Ohms")

# --- Step 4: Apply the startup margin to find the target impedance ---
# The target impedance seen by the diode (RL_prime) includes the margin
RL_prime = Rs_mag * (1 + margin)

print(f"\nFor optimal power transfer with a {margin*100:.0f}% startup margin, the diode must see a load (RL') of:")
print(f"RL' = |Rs| * (1 + margin) = {Rs_mag:.4f} Ohms * {1+margin} = {RL_prime:.4f} Ohms")

# --- Step 5: Calculate the final impedance transformation ratio ---
# The ratio transforms the 50 ohm load to RL_prime
# Ratio = RL_load / RL_prime
transformation_ratio = RL_load / RL_prime

print(f"\nThe impedance transformation ratio from the {RL_load} Ohm load to the diode is:")
print(f"Ratio = Load Impedance / Diode-side Impedance = {RL_load} Ohms / {RL_prime:.4f} Ohms")
print(f"Final Transformation Ratio: {transformation_ratio:.4f}")

print(f"<<<{transformation_ratio}>>>")