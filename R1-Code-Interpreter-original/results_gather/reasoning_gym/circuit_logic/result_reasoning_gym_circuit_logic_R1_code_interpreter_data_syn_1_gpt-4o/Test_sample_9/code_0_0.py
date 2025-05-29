# Given inputs
A = 1
B = 1
C = 0
D = 0
E = 0
F = 1

# Negate operations
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F

# Evaluate the circuit step by step
# AND gates
and1 = A and B
and2 = C and D
and3 = E and F

# NAND gates
nand1 = 1 - (A and B)
nand2 = 1 - (C and D)
nand3 = 1 - (E and F)

# XOR gates
xor1 = A ^ B
xor2 = C ^ D
xor3 = E ^ F

# Combine the results based on the circuit
# The final output is determined by the last XOR gate
final_output = xor1 ^ xor2 ^ xor3

# Print the final output
print(final_output)