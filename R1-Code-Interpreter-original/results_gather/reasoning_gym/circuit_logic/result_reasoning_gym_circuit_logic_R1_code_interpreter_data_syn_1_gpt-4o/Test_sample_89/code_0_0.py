# Initial inputs
A = 0
B = 0
C = 0

# XOR gate evaluations
XOR1 = A ^ B
XOR2 = B ^ C
XOR3 = C ^ A

# Negate outputs
NOT_XOR1 = not XOR1
NOT_XOR2 = not XOR2
NOT_XOR3 = not XOR3

# XOR gate evaluations after negation
XOR4 = NOT_XOR1 ^ NOT_XOR2
XOR5 = NOT_XOR2 ^ NOT_XOR3

# AND gate evaluations
AND1 = A and B
AND2 = NOT_XOR1 and NOT_XOR2

# NAND gate evaluations
NAND1 = not (XOR4 and XOR5)

# NOR gate evaluations
NOR1 = not (XOR5 or NAND1)

# Final AND gate
OUT = AND1 and AND2 and NOR1

print(OUT)