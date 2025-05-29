# Input values
A, B, C, D, E, F, G = 1, 1, 1, 1, 1, 0, 0

# Negations
not_G = 1 - G
not_F = 1 - F
not_E = 1 - E
not_D = 1 - D
not_C = 1 - C
not_B = 1 - B

# AND gates
AND1 = not_G & not_F
AND2 = not_E & not_D
AND3 = not_C & not_B

# NAND gates
NAND1 = 1 - (not_G & not_F)
NAND2 = 1 - (not_E & not_D)
NAND3 = 1 - (not_C & not_B)

# XOR gates
XOR1 = A ^ AND1
XOR2 = B ^ NAND1
XOR3 = C ^ AND2
XOR4 = D ^ NAND2
XOR5 = E ^ AND3
XOR6 = F ^ NAND3

# Final XOR for OUT
OUT = XOR1 ^ XOR2 ^ XOR3 ^ XOR4 ^ XOR5 ^ XOR6

print(OUT)