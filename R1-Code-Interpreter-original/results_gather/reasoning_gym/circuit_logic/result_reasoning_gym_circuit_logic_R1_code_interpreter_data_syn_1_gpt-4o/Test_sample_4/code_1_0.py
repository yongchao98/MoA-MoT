# Input values
A = 1
B = 0
C = 0
D = 1
E = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E

# AND gates
AND1 = A and B
AND2 = C and D
AND3 = not_A and not_B
AND4 = not_C and not_D
AND5 = not_E and E

# XOR gates
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = not_A ^ not_B
XOR4 = not_C ^ not_D
XOR5 = not_E ^ E

# NAND gates
NAND1 = 1 - (A and B)
NAND2 = 1 - (C and D)
NAND3 = 1 - (not_A and not_B)
NAND4 = 1 - (not_C and not_D)
NAND5 = 1 - (not_E and E)

# OR gates
OR1 = A or B
OR2 = C or D
OR3 = not_A or not_B
OR4 = not_C or not_D
OR5 = not_E or E

# Final output calculation
OUT = (AND1 or AND2 or AND3 or AND4 or AND5) and (XOR1 or XOR2 or XOR3 or XOR4 or XOR5) and (NAND1 or NAND2 or NAND3 or NAND4 or NAND5) and (OR1 or OR2 or OR3 or OR4 or OR5)

print(OUT)