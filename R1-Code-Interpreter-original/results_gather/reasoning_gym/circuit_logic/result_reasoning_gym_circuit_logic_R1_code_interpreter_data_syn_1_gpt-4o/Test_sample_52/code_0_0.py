# Input assignments
A = 0
B = 0
C = 0
D = 1
E = 0
F = 1
G = 1
H = 1
I = 0
J = 1

# Negation (NOT) gate function
def NOT(x):
    return 1 if x == 0 else 0

# AND gate function
def AND(*args):
    return 1 if all(args) else 0

# NAND gate function
def NAND(*args):
    return NOT(AND(*args))

# XOR gate function
def XOR(*args):
    return 1 if sum(args) % 2 == 1 else 0

# Evaluate the circuit step by step
# Following the connections in the circuit diagram

# First layer of NOT gates
not_J = NOT(J)
not_I = NOT(I)
not_H = NOT(H)
not_G = NOT(G)
not_F = NOT(F)
not_E = NOT(E)
not_D = NOT(D)
not_C = NOT(C)
not_B = NOT(B)
not_A = NOT(A)

# Second layer of AND gates
and_1 = AND(not_J, not_I)
and_2 = AND(not_H, not_G)
and_3 = AND(not_F, not_E)
and_4 = AND(not_D, not_C)
and_5 = AND(not_B, not_A)

# Third layer of XOR gates
xor_1 = XOR(and_1, and_2)
xor_2 = XOR(and_3, and_4)

# Fourth layer of NAND gate
nand_1 = NAND(xor_1, xor_2)

# Fifth layer of AND gate
and_6 = AND(and_5, nand_1)

# Final XOR gate for output
output = XOR(and_6, xor_1)

# Print the final output
print(output)