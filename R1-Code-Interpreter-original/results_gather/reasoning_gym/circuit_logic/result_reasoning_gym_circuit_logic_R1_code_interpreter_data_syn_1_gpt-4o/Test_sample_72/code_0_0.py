# Input values
A = 0
B = 1
C = 1

# XOR gates
XOR1 = A ^ B
XOR2 = B ^ C

# Negation of XOR2
NOT_XOR2 = 1 if XOR2 == 0 else 0

# AND gates
AND1 = XOR1 & NOT_XOR2
AND2 = XOR2 & XOR1

# NOR gates
NOR1 = 1 if (AND1 == 0 and AND2 == 0) else 0
NOR2 = 1 if (NOR1 == 0 and XOR2 == 0) else 0

# Final output
OUT = NOR2

print(OUT)