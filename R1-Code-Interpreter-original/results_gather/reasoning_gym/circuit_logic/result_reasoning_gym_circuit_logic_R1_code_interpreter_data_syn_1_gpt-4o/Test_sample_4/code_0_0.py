# Input values
A = 1
B = 0
C = 0
D = 1
E = 1

# Negations
not_A = not A
not_B = not B
not_C = not C
not_D = not D
not_E = not E

# AND Gates
AND1 = A and B
AND2 = C and D
AND3 = E and not_A
AND4 = not_B and not_C
AND5 = not_D and not_E

# XOR Gates
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = E ^ not_A
XOR4 = not_B ^ not_C
XOR5 = not_D ^ not_E

# NAND Gates
NAND1 = not (A and B)
NAND2 = not (C and D)
NAND3 = not (E and not_A)
NAND4 = not (not_B and not_C)
NAND5 = not (not_D and not_E)

# OR Gates
OR1 = AND1 or AND2
OR2 = AND3 or AND4
OR3 = AND5 or XOR1
OR4 = XOR2 or XOR3
OR5 = XOR4 or XOR5
OR6 = NAND1 or NAND2
OR7 = NAND3 or NAND4
OR8 = NAND5 or OR1

# Final Output
final_output = OR2 or OR3 or OR4 or OR5 or OR6 or OR7 or OR8

print(final_output)