# Input values
A = 1
B = 1
C = 0
D = 0
E = 1

# Negate
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# AND gates
AND1 = A and B
AND2 = C and D
AND3 = E and not_A
AND4 = not_B and not_C
AND5 = not_D and not_E

# NAND gates
NAND1 = 1 - (A and B)
NAND2 = 1 - (C and D)
NAND3 = 1 - (E and not_A)
NAND4 = 1 - (not_B and not_C)
NAND5 = 1 - (not_D and not_E)

# XOR gates
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = E ^ not_A
XOR4 = not_B ^ not_C
XOR5 = not_D ^ not_E

# OR gates
OR1 = AND1 or AND2
OR2 = AND3 or AND4
OR3 = AND5 or NAND1
OR4 = NAND2 or NAND3
OR5 = NAND4 or NAND5
OR6 = XOR1 or XOR2
OR7 = XOR3 or XOR4
OR8 = XOR5 or OR1

# Final output
OUT = OR2 or OR3 or OR4 or OR5 or OR6 or OR7 or OR8

print(OUT)