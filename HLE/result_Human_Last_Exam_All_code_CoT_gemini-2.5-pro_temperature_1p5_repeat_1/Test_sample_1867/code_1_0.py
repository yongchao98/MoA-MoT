import math

# Define given constants
Io = 1e-9  # Reverse saturation current in Amperes
n = 1.5    # Diode ideality factor
T = 300    # Ambient temperature in Kelvin
V1 = 0.78  # Start voltage of linear region in Volts
V2 = 0.98  # End voltage of linear region in Volts
I2 = 0.445 # Current at V2 in Amperes
R_L = 50   # Load resistance in Ohms
margin = 0.20 # Startup margin (20%)

# Physical constants
k = 1.380649e-23  # Boltzmann constant in J/K
q = 1.60217663e-19  # Electron charge in Coulombs

# Step 1: Calculate the thermal voltage (Vt)
Vt = (k * T) / q
print(f"Step 1: Calculated Thermal Voltage (Vt)")
print(f"Vt = (k * T) / q = ({k:.6e} * {T}) / {q:.6e} = {Vt:.5f} V\n")

# Step 2: Calculate the current I1 at V1 using the Shockley diode equation
I1 = Io * (math.exp(V1 / (n * Vt)) - 1)
print(f"Step 2: Calculated Diode Current at V1 (I1)")
print(f"I1 = Io * (exp(V1 / (n * Vt)) - 1)")
print(f"I1 = {Io:.1e} * (exp({V1} / ({n} * {Vt:.5f})) - 1) = {I1:.5f} A\n")

# Step 3: Calculate the dynamic resistance (rd) in the linear region
# rd = dV/dI = (V2 - V1) / (I2 - I1)
rd = (V2 - V1) / (I2 - I1)
print(f"Step 3: Calculated Diode Dynamic Resistance (rd)")
print(f"rd = (V2 - V1) / (I2 - I1)")
print(f"rd = ({V2} - {V1}) / ({I2} - {I1:.5f}) = {rd:.5f} Ohms\n")

# Step 4: Determine the target load impedance for the diode with a startup margin
# For oscillation startup, R_transformed < |rd|.
# With a margin, R_transformed = |rd| / (1 + margin)
abs_rd = abs(rd)
R_transformed = abs_rd / (1 + margin)
print(f"Step 4: Calculated Target Load Impedance (R_transformed)")
print(f"R_transformed = |rd| / (1 + margin)")
print(f"R_transformed = {abs_rd:.5f} / (1 + {margin}) = {R_transformed:.5f} Ohms\n")

# Step 5: Calculate the impedance transformation ratio (K)
# R_transformed = K * R_L
K = R_transformed / R_L
print(f"Step 5: Calculated Impedance Transformation Ratio (K)")
print(f"K = R_transformed / R_L")
print(f"K = {R_transformed:.5f} / {R_L} = {K:.5f}\n")

# Final Answer
print("The final impedance transformation ratio should be:")
print(f"<<<{K:.5f}>>>")