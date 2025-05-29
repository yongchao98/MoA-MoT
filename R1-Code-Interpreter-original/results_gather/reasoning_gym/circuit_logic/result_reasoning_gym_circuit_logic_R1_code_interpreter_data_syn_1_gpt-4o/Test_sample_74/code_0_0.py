# Input values
A = 0
B = 1
C = 1
D = 1
E = 0
F = 1

# Negate gates
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F

# AND gates
AND1 = not_A and not_B
AND2 = not_C and not_D
AND3 = not_E and not_F
AND4 = A and B
AND5 = C and D
AND6 = E and F

# NAND gates
NAND1 = 1 - (not_A and not_B)
NAND2 = 1 - (not_C and not_D)
NAND3 = 1 - (not_E and not_F)
NAND4 = 1 - (A and B)
NAND5 = 1 - (C and D)
NAND6 = 1 - (E and F)

# XOR gates
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = E ^ F

# Final AND gates
AND7 = NAND1 and NAND2
AND8 = NAND3 and NAND4
AND9 = NAND5 and NAND6
AND10 = XOR1 and XOR2
AND11 = XOR3 and AND7
AND12 = AND8 and AND9
AND13 = AND10 and AND11
AND14 = AND12 and AND13

# Output
print(AND14)