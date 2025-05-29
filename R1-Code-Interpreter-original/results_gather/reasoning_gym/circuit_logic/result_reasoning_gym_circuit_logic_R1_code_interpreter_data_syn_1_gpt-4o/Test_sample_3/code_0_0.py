# Input values
A = 0
B = 1
C = 1
D = 0
E = 0

# AND gates
AND1 = D and E
AND2 = C and AND1
AND3 = B and AND2
AND4 = A and AND3

# NAND gates
NAND1 = not (D and E)
NAND2 = not (C and NAND1)
NAND3 = not (B and NAND2)
NAND4 = not (A and NAND3)

# NOT gates
NOT1 = not NAND1
NOT2 = not NAND2
NOT3 = not NAND3
NOT4 = not NAND4

# XOR gates
XOR1 = NOT1 != NOT2
XOR2 = NOT3 != NOT4
OUT = XOR1 != XOR2

print(OUT)