# Input values
A = 0
B = 0
C = 0
D = 0
E = 1
F = 1
G = 0
H = 0

# Negate gates
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G
not_H = 1 - H

# NAND gates
nand_1 = 1 - (H * G)
nand_2 = 1 - (F * E)
nand_3 = 1 - (D * C)
nand_4 = 1 - (B * A)

# XOR gates
xor_1 = (nand_1 + nand_2) % 2
xor_2 = (nand_3 + nand_4) % 2

# AND gates
and_1 = nand_1 * nand_2
and_2 = nand_3 * nand_4
and_3 = xor_1 * xor_2

# Final AND gate for output
OUT = and_1 * and_2 * and_3

print(OUT)