# Input values
A = 1
B = 1
C = 1
D = 0
E = 1
F = 1
G = 0

# Negate operations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G

# AND operations
and1 = A and B
and2 = C and D
and3 = E and F
and4 = G and not_A

# NAND operations
nand1 = 1 - (A and B)
nand2 = 1 - (C and D)
nand3 = 1 - (E and F)
nand4 = 1 - (G and not_A)

# XOR operations
xor1 = A ^ B
xor2 = C ^ D
xor3 = E ^ F
xor4 = G ^ not_A

# Final output calculation based on the circuit
# The final output is determined by the specific connections in the circuit
# Assuming the final output is a combination of the above results
# This needs to be adjusted based on the specific circuit logic

# Example final output calculation (this is a placeholder, adjust as needed)
final_output = xor1 ^ xor2 ^ xor3 ^ xor4 ^ nand1 ^ nand2 ^ nand3 ^ nand4

print(final_output)