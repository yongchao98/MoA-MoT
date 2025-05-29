# Define the logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def XOR(x, y):
    return x ^ y

def NOT(x):
    return ~x & 1

def NOR(x, y):
    return ~(x | y) & 1

# Given inputs
A = 0
B = 0
C = 1
D = 1
E = 1

# Evaluate the circuit step by step
# Following the circuit diagram from bottom to top

# First layer of gates
negate_E = NOT(E)
nand_D_E = NAND(D, E)
nand_C_D = NAND(C, D)
nand_B_C = NAND(B, C)
nand_A_B = NAND(A, B)

# Second layer of gates
xor_negate_E_nand_D_E = XOR(negate_E, nand_D_E)
xor_nand_C_D_nand_B_C = XOR(nand_C_D, nand_B_C)
xor_nand_A_B_negate_E = XOR(nand_A_B, negate_E)

# Third layer of gates
nand_xor_negate_E_nand_D_E_xor_nand_C_D_nand_B_C = NAND(xor_negate_E_nand_D_E, xor_nand_C_D_nand_B_C)
nand_xor_nand_A_B_negate_E_nand_xor_negate_E_nand_D_E_xor_nand_C_D_nand_B_C = NAND(xor_nand_A_B_negate_E, nand_xor_negate_E_nand_D_E_xor_nand_C_D_nand_B_C)

# Fourth layer of gates
nor_nand_xor_nand_A_B_negate_E_nand_xor_negate_E_nand_D_E_xor_nand_C_D_nand_B_C_nand_xor_negate_E_nand_D_E = NOR(nand_xor_nand_A_B_negate_E_nand_xor_negate_E_nand_D_E_xor_nand_C_D_nand_B_C, xor_negate_E_nand_D_E)

# Final output
OUT = nor_nand_xor_nand_A_B_negate_E_nand_xor_negate_E_nand_D_E_xor_nand_C_D_nand_B_C_nand_xor_negate_E_nand_D_E

print(OUT)