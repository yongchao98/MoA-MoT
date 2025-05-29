# Input values
A, B, C, D, E, F, G = 0, 1, 0, 1, 0, 0, 0

# Negate the inputs
A_neg = 1 - A
B_neg = 1 - B
C_neg = 1 - C
D_neg = 1 - D
E_neg = 1 - E
F_neg = 1 - F
G_neg = 1 - G

# Evaluate the circuit
# Following the circuit diagram, evaluate step by step
nand1 = 1 - (A_neg & B_neg)
nand2 = 1 - (C_neg & D_neg)
nand3 = 1 - (E_neg & F_neg)
nand4 = 1 - (nand1 & nand2)
nand5 = 1 - (nand3 & G_neg)

xor1 = nand4 ^ nand5
xor2 = xor1 ^ nand2
xor3 = xor2 ^ nand3

# Final output
output = xor3
print(output)