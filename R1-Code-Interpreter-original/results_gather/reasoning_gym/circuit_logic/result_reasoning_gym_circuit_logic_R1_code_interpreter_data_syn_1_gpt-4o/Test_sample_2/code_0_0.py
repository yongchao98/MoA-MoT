# Input values
A = 0
B = 0
C = 1
D = 1
E = 1
F = 0

# Negate gates
negate_F = 1 if F == 0 else 0
negate_D = 1 if D == 0 else 0

# AND gates
AND1 = C and D
AND2 = B and AND1
AND3 = negate_D and E

# NAND gates
NAND1 = 1 if not (E and F) else 0
NAND2 = 1 if not (negate_F and NAND1) else 0

# XOR gates
XOR1 = A ^ B
XOR2 = XOR1 ^ C
XOR3 = XOR2 ^ AND2
XOR4 = XOR3 ^ NAND2
XOR5 = XOR4 ^ AND3

# Final XOR
final_output = XOR5 ^ NAND1

print(final_output)