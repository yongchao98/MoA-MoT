import math

# Description:
# This script calculates the required impedance transformation ratio for a diode
# acting as a signal source connected to a 50 ohm load.

# --- Given Parameters ---
Io = 1e-9       # Reverse saturation current in Amperes
n = 1.5         # Diode ideality factor
T = 300         # Ambient temperature in Kelvin
V1 = 0.78       # Start voltage of linear region in Volts
V2 = 0.98       # End voltage of linear region in Volts
I2 = 0.445      # Current at V2 in Amperes
RL = 50.0       # Load resistance in Ohms
margin = 0.20   # Startup margin (20%)

# --- Physical Constants ---
k = 1.380649e-23  # Boltzmann's constant in J/K
q = 1.602176634e-19  # Elementary charge in Coulombs

# Step 1: Calculate the thermal voltage (Vt)
Vt = (k * T) / q

# Step 2: Calculate the current I1 at the start of the linear region (V1)
# We assume the linear behavior starts at the point defined by the standard diode equation.
exponent = V1 / (n * Vt)
I1 = Io * (math.exp(exponent) - 1)

# Step 3: Calculate the dynamic resistance (rd) of the diode in the linear region.
# A negative value indicates the diode is behaving as an active signal source.
delta_V = V2 - V1
delta_I = I2 - I1
rd = delta_V / delta_I

# Step 4: Calculate the optimal impedance transformation ratio for maximum power transfer.
# The ratio is |rd| / RL, which scales the load impedance to match the source.
optimal_ratio = abs(rd) / RL

# Step 5: Apply the 20% startup margin
final_ratio = optimal_ratio * (1 + margin)

# --- Output the results ---
print("Calculating the impedance transformation ratio...")
print("-" * 50)
print(f"Calculated dynamic source resistance (rd): {rd:.3f} Ohms")
print(f"Load resistance (RL): {RL:.1f} Ohms")
print(f"Optimal impedance transformation ratio (|rd| / RL): {optimal_ratio:.4f}")
print(f"Required startup margin: {margin*100:.0f}%")
print("-" * 50)
print("The final transformation ratio is calculated as follows:")
print(f"Final Ratio = (Optimal Ratio) * (1 + Margin)")
print(f"Final Ratio = {optimal_ratio:.4f} * (1 + {margin})")
print(f"Final Ratio = {final_ratio:.4f}")
print("-" * 50)

print(f"<<<{final_ratio:.4f}>>>")