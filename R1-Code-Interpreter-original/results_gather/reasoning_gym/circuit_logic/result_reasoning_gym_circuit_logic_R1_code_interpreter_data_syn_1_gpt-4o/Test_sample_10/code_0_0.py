# Input values
A = 1
B = 0
C = 1
D = 1
E = 0
F = 1
G = 0
H = 0
I = 0

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G
not_H = 1 - H
not_I = 1 - I

# AND gates
AND1 = not_A and not_B
AND2 = not_C and not_D
AND3 = not_E and not_F
AND4 = not_G and not_H
AND5 = not_I and A
AND6 = B and C
AND7 = D and E
AND8 = F and G
AND9 = H and I

# NAND gates
NAND1 = 1 - (AND4 and AND5)
NAND2 = 1 - (AND6 and AND7)
NAND3 = 1 - (AND8 and AND9)

# XOR gates
XOR1 = AND1 ^ AND2
XOR2 = AND3 ^ NAND1
XOR3 = NAND2 ^ NAND3

# OR gates
OR1 = XOR1 or XOR2
OR2 = OR1 or XOR3

# Output
print(OR2)