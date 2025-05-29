# Input values
A = 0
B = 1
C = 1
D = 0
E = 0
F = 0

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F

# XOR Gates
XOR1 = not_F ^ not_E
XOR2 = not_D ^ not_C
XOR3 = not_B ^ not_A
XOR4 = D ^ C
XOR5 = E ^ F

# AND Gates
AND1 = XOR1 & XOR2
AND2 = XOR3 & XOR4

# NAND Gates
NAND1 = 1 - (XOR5 & B)
NAND2 = 1 - (A & XOR5)

# OR Gates
OR1 = AND1 | NAND1
OR2 = AND2 | NAND2
OR3 = OR1 | OR2

# Final Output
print(OR3)