# Input values
A = 0
B = 0
C = 0
D = 0
E = 1
F = 0
G = 0

# Negate values
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G

# XOR gates
XOR1 = A ^ B
XOR2 = C ^ D
XOR3 = E ^ F
XOR4 = G ^ not_G

# AND gates
AND1 = XOR1 & XOR2
AND2 = XOR3 & XOR4
AND3 = AND1 & AND2

# OR gates
OR1 = AND3 | not_A
OR2 = OR1 | not_B
OR3 = OR2 | not_C
OR4 = OR3 | not_D
OR5 = OR4 | not_E
OR6 = OR5 | not_F
OR7 = OR6 | not_G

# Final output
print(OR7)