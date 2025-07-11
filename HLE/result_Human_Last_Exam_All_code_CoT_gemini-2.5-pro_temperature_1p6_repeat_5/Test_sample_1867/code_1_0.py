import math

# --- Step 1: Define Constants and Initial Values ---
# Boltzman's Constant (J/K)
k = 1.380649e-23
# Elementary Charge (C)
q = 1.602176634e-19
# Temperature (K)
T = 300.0
# Reverse Saturation Current (A)
Io = 1e-9
# Diode Ideality Factor
n = 1.5
# Diode voltages for the linear region (V)
V1 = 0.78
V2 = 0.98
# Diode current at V2 (A)
I2 = 0.445
# Load Resistance (ohms)
RL = 50.0
# Startup margin
margin = 0.20

# --- Step 2: Calculate Thermal Voltage (Vt) ---
Vt = (k * T) / q
print(f"Step 1: Calculating Thermal Voltage (Vt)")
print(f"Vt = (k * T) / q = ({k:.4e} * {T}) / {q:.4e} = {Vt:.5f} V\n")

# --- Step 3: Calculate Diode Current I1 at V1 ---
exponent = V1 / (n * Vt)
# The -1 is negligible at this forward bias but included for completeness
I1 = Io * (math.exp(exponent) - 1)
print(f"Step 2: Calculating Diode Current (I1) at V1 = {V1} V")
print(f"I1 = Io * (exp(V1 / (n * Vt)) - 1)")
print(f"I1 = {Io:.1e} * (exp({V1} / ({n} * {Vt:.5f})) - 1) = {I1:.5f} A\n")

# --- Step 4: Calculate Dynamic Source Resistance (Rs) ---
delta_V = V2 - V1
delta_I = I2 - I1
Rs = delta_V / delta_I
print(f"Step 3: Calculating Dynamic Source Resistance (Rs)")
print(f"Rs = (V2 - V1) / (I2 - I1) = ({V2} - {V1}) / ({I2:.5f} - {I1:.5f})")
print(f"Rs = {delta_V:.2f} V / {delta_I:.5f} A = {Rs:.5f} ohms\n")

# --- Step 5: Determine Optimal Load Impedance for Maximum Power Transfer ---
# For a negative resistance source, the optimal load is -Rs
RL_prime_optimal = -Rs
print(f"Step 4: Determining Optimal Load Impedance (RL') for the Diode")
print(f"For max power transfer from a negative resistance source, RL'_optimal = -Rs")
print(f"RL'_optimal = -({Rs:.5f}) = {RL_prime_optimal:.5f} ohms\n")

# --- Step 6: Apply the 20% Startup Margin ---
RL_prime_adjusted = RL_prime_optimal * (1 - margin)
print(f"Step 5: Applying 20% Startup Margin")
print(f"For reliable startup, the target load must be less than the optimal value.")
print(f"RL'_adjusted = RL'_optimal * (1 - margin) = {RL_prime_optimal:.5f} * (1 - {margin}) = {RL_prime_adjusted:.5f} ohms\n")

# --- Step 7: Calculate the Final Impedance Transformation Ratio ---
ratio = RL_prime_adjusted / RL
print(f"Step 6: Calculating the Impedance Transformation Ratio")
print(f"The transformer must make the {RL} ohm load appear as {RL_prime_adjusted:.5f} ohms.")
print(f"Ratio = RL'_adjusted / RL = {RL_prime_adjusted:.5f} / {RL} = {ratio:.5f}\n")

print("The final impedance transformation ratio should be:")
print(f"<<<{ratio:.5f}>>>")