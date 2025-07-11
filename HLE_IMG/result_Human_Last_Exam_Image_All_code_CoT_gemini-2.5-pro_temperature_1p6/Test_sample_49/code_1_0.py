import math

# This script prints the formula for the angular cutoff frequency at node a0
# in the given ladder network.

# The derivation involves finding the Thevenin equivalent resistance (R_th)
# seen by the capacitor C. The analysis of the infinite ladder leads to:
# R_th = (1 + sqrt(3)) * r

# The angular cutoff frequency is then given by omega_c = 1 / (R_th * C).

# The final formula expresses the frequency in terms of the circuit's
# resistance 'r' and capacitance 'C'. The numbers in the equation are 1 and 3.
num_one = 1
num_three = 3

print("The angular cutoff frequency omega_c at node a0 is given by the equation:")
print(f"omega_c = {num_one} / (({num_one} + sqrt({num_three})) * r * C)")