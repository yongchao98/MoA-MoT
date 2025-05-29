# Define the input values
A = 0
B = 0
C = 1
D = 0
E = 1
F = 1
G = 1

# Negate operations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G

# NAND operations
nand1 = 1 - (A & B)
nand2 = 1 - (C & D)
nand3 = 1 - (E & F)
nand4 = 1 - (G & nand3)

# XOR operations
xor1 = (nand1 ^ not_G)
xor2 = (nand2 ^ not_F)
xor3 = (nand4 ^ not_E)
xor4 = (not_D ^ not_C)
xor5 = (not_B ^ not_A)

# AND operations
and1 = xor1 & xor2
and2 = xor3 & xor4
and3 = and1 & and2
and4 = and3 & xor5

# Final output
print(and4)