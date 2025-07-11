import math

# Step 1: Define constants from the problem statement and physics
Io = 1e-9      # Reverse saturation current in Amperes
n = 1.5        # Diode ideality factor
T = 300        # Ambient temperature in Kelvin
V1 = 0.78      # Start voltage of the linear region in Volts
V2 = 0.98      # End voltage of the linear region in Volts
I2 = 0.445     # Current at V2 in Amperes
RL = 50.0      # Load resistance in Ohms
margin = 0.20  # Startup margin

# Physical constants
k = 1.380649e-23  # Boltzmann's constant in J/K
q = 1.60217663e-19  # Elementary charge in C

# Step 2: Calculate the thermal voltage (Vt)
Vt = (k * T) / q

# Step 3: Calculate the current I1 at voltage V1 using the diode equation
# The '-1' term is negligible for this forward bias but included for accuracy.
I1 = Io * (math.exp(V1 / (n * Vt)) - 1)

# Step 4: Calculate the dynamic resistance (rd) in the linear region
# rd = dV / dI = (V2 - V1) / (I2 - I1)
delta_V = V2 - V1
delta_I = I2 - I1
rd = delta_V / delta_I

# Step 5: Determine the target transformed load resistance (RL_prime)
# For optimum power transfer, the transformed load resistance must match the
# magnitude of the source's dynamic resistance.
RL_prime = abs(rd)

# Step 6: Calculate the ideal impedance transformation ratio (N_ideal)
# N = Z_source_side / Z_load_side = RL_prime / RL
N_ideal = RL_prime / RL

# Step 7: Apply the 20% startup margin
N_final = N_ideal * (1 + margin)

# Step 8: Print the final calculated values to show the result
# The "final equation" is the series of calculations performed.
# Here are the values that constitute the final result:
print(f"Diode dynamic resistance (rd): {rd:.5f} Ohms")
print(f"Required transformed load resistance (RL'): {RL_prime:.5f} Ohms")
print(f"Ideal transformation ratio: {N_ideal:.5f}")
print(f"Final transformation ratio with 20% margin: {N_final:.5f}")
