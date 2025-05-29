# Input assignments
A = 1
B = 0
C = 1
D = 0
E = 0
F = 1
G = 0
H = 1

# Negation function
def negate(x):
    return 1 if x == 0 else 0

# NAND function
def nand(x, y):
    return negate(x & y)

# XOR function
def xor(x, y):
    return x ^ y

# Step 1: Evaluate the negations and NAND gates
neg_H = negate(H)
nand_GH = nand(G, H)
nand_FG = nand(F, G)
nand_EF = nand(E, F)
nand_DE = nand(D, E)
nand_CD = nand(C, D)
nand_BC = nand(B, C)
nand_AB = nand(A, B)

# Step 2: Evaluate the XOR gates
xor_1 = xor(nand_GH, nand_FG)
xor_2 = xor(nand_EF, nand_DE)
xor_3 = xor(nand_CD, nand_BC)
xor_4 = xor(nand_AB, neg_H)

# Step 3: Continue evaluating the circuit
xor_5 = xor(xor_1, xor_2)
xor_6 = xor(xor_3, xor_4)
final_output = xor(xor_5, xor_6)

# Print the final output
print(final_output)