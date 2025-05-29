# Input values
A = 1
B = 0
C = 1
D = 1
E = 0
F = 1
G = 1
H = 1

# Negate gates
not_A = not A
not_B = not B
not_C = not C
not_D = not D
not_E = not E
not_F = not F
not_G = not G
not_H = not H

# NAND gates
NAND1 = not (not_H and G)
NAND2 = not (NAND1 and F)
NAND3 = not (NAND2 and E)
NAND4 = not (NAND3 and D)
NAND5 = not (NAND4 and C)
NAND6 = not (NAND5 and B)
NAND7 = not (NAND6 and A)

# XOR gates
XOR1 = NAND7 ^ NAND6
XOR2 = XOR1 ^ NAND5
XOR3 = XOR2 ^ NAND4
XOR4 = XOR3 ^ NAND3

# AND gates
AND1 = XOR4 and NAND2
AND2 = AND1 and NAND1

# Final output
OUT = AND2

print(OUT)