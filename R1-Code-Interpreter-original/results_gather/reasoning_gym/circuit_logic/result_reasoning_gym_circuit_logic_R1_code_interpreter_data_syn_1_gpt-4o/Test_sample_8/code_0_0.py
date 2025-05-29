# Input values
A = 1
B = 1
C = 1
D = 0
E = 0
F = 1
G = 1
H = 0

# Negations
H_neg = not H
G_neg = not G
F_neg = not F
E_neg = not E
D_neg = not D
C_neg = not C
B_neg = not B
A_neg = not A

# AND gates
AND1 = H_neg and G
AND2 = F_neg and E
AND3 = D_neg and C
AND4 = B_neg and A

# XOR gates
XOR1 = AND1 ^ AND2
XOR2 = AND3 ^ AND4
XOR3 = XOR1 ^ XOR2

# OR gates
OR1 = XOR3 or AND2
OR2 = OR1 or AND4

# Final output
output = OR2
print(output)