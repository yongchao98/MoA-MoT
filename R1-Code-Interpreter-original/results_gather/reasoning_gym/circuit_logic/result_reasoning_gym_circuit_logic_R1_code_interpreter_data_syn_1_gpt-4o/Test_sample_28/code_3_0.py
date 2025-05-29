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

# First layer of negations
neg1 = 1 - A
neg2 = 1 - B
neg3 = 1 - C
neg4 = 1 - D
neg5 = 1 - E
neg6 = 1 - F
neg7 = 1 - G

# NAND gates
nand1 = 1 - (neg1 & neg2)
nand2 = 1 - (neg3 & neg4)
nand3 = 1 - (neg5 & neg6)
nand4 = 1 - (nand1 & nand2)
nand5 = 1 - (nand3 & neg7)

# XOR gates
xor1 = nand4 ^ nand5
xor2 = xor1 ^ nand2
xor3 = xor2 ^ nand3

# Final output
output = xor3
print(output)