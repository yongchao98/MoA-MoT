# Input values
A = 0
B = 0
C = 1
D = 0

# Negate gates
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D

# NAND gates
nand1 = 1 - (not_A & not_B)
nand2 = 1 - (not_C & not_D)
nand3 = 1 - (A & B)
nand4 = 1 - (C & D)

# XOR gates
xor1 = nand1 ^ nand2
xor2 = nand3 ^ nand4

# AND gate
out = xor1 & xor2

print(out)