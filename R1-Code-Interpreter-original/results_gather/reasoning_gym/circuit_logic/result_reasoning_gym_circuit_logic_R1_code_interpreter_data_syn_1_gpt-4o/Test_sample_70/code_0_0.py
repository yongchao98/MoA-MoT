# Input values
A = 1
B = 0
C = 1
D = 1
E = 1
F = 0
G = 0
H = 1

# Negation
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# AND gates
and1 = H and G
and2 = F and and1
and3 = E and and2
and4 = D and and3

# XOR gates
xor1 = B ^ not_C
xor2 = xor1 ^ not_D
xor3 = xor2 ^ not_E
xor4 = A ^ xor3

# NAND gates
nand1 = 1 - (and4 and xor4)
nand2 = 1 - (nand1 and xor4)

# Final XOR gate
final_output = nand2 ^ xor4

print(final_output)