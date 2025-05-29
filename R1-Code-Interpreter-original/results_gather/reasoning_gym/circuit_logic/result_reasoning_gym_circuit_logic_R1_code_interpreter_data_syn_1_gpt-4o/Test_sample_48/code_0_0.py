# Input values
A = 1
B = 0
C = 1
D = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D

# XOR gates
XOR1 = not_D ^ C
XOR2 = not_C ^ B
XOR3 = not_B ^ A
XOR4 = not_A ^ D

# AND gates
AND1 = XOR1 & XOR2
AND2 = XOR3 & XOR4
AND3 = not_C & not_D
AND4 = not_B & not_A

# NAND gate
NAND1 = 1 - (AND3 & AND4)

# Final AND gates
Final_AND1 = AND1 & AND2
Final_AND2 = Final_AND1 & NAND1

# Output
print(Final_AND2)