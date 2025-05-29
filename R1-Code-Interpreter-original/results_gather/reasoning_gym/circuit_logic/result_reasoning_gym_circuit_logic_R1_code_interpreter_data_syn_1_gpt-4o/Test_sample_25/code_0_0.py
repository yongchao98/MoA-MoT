# Input values
A = 1
B = 0
C = 1
D = 1
E = 0
F = 0
G = 1
H = 1

# Negations
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D
not_E = 1 - E
not_F = 1 - F
not_G = 1 - G
not_H = 1 - H

# XOR operations
XOR1 = not_H ^ G
XOR2 = not_G ^ F
XOR3 = not_F ^ E
XOR4 = not_E ^ D
XOR5 = not_D ^ C
XOR6 = not_C ^ B
XOR7 = not_B ^ A

# Further XOR operations
XOR8 = XOR1 ^ XOR2
XOR9 = XOR3 ^ XOR4
XOR10 = XOR5 ^ XOR6
XOR11 = XOR7 ^ not_A

# Final XOR operations
XOR12 = XOR8 ^ XOR9
XOR13 = XOR10 ^ XOR11

# Final output
OUT = XOR12 ^ XOR13

print(OUT)