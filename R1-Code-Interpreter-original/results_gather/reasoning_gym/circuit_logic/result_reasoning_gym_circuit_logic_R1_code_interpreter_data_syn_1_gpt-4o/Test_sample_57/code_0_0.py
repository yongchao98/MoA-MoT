# Input values
A = 0
B = 0
C = 0

# Negation
not_C = 1 if C == 0 else 0

# AND gates
AND1 = A and B
AND2 = B and C
AND3 = A and C
AND4 = A and not_C
AND5 = B and (1 if AND1 == 0 else 0)

# XOR gates
XOR1 = AND2 ^ AND3
XOR2 = AND4 ^ AND5
OUT = XOR1 ^ XOR2

print(OUT)