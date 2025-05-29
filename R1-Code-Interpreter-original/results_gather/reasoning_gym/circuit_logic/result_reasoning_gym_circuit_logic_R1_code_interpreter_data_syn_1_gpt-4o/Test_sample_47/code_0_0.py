# Input values
A = 1
B = 1
C = 1
D = 0
E = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# NAND gates
NAND1 = 1 - (not_A and not_B)
NAND2 = 1 - (not_C and not_D)
NAND3 = 1 - (not_E and not_A)
NAND4 = 1 - (not_B and not_C)

# AND gates
AND1 = NAND1 and NAND2
AND2 = NAND3 and NAND4

# XOR gates
XOR1 = not_A ^ not_B
XOR2 = not_C ^ not_D
XOR3 = XOR1 ^ XOR2

# Final output
OUT = AND1 and AND2 and XOR3

print(OUT)