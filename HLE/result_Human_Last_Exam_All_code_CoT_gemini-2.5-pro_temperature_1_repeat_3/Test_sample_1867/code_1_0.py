import math

# Step 1: Define the given parameters from the problem description.
v1 = 0.78  # Start voltage in Volts
v2 = 0.98  # End voltage in Volts
delta_i = 0.445  # The change in current over the voltage range, in Amps.
rl = 50.0  # Load resistance in Ohms
margin = 0.20  # 20% margin

# Step 2: Calculate the change in voltage (delta_v).
delta_v = v2 - v1

# Step 3: Calculate the dynamic resistance of the diode (the source resistance, Rs).
# This is calculated as the change in voltage divided by the change in current.
rs = delta_v / delta_i

# Step 4: Calculate the ideal impedance transformation ratio for maximum power transfer.
# This ratio should make the load impedance appear equal to the source impedance.
# Ratio = Z_source / Z_load = Rs / RL
ideal_ratio = rs / rl

# Step 5: Apply the 20% startup margin to the ideal ratio.
final_ratio = ideal_ratio * (1 + margin)

# Step 6: Print the results, showing the values used in the final calculation.
print(f"The change in voltage (ΔV) is: {v2} V - {v1} V = {delta_v:.2f} V")
print(f"The change in current (ΔI) is given as: {delta_i} A")
print(f"The diode's dynamic resistance (Rs = ΔV/ΔI) is: {rs:.4f} Ohms")
print(f"The load resistance (RL) is: {rl} Ohms")
print(f"The ideal impedance transformation ratio (Rs/RL) is: {ideal_ratio:.5f}")
print(f"The final required impedance transformation ratio with a 20% margin is ({ideal_ratio:.5f} * {1 + margin}) = {final_ratio:.5f}")
