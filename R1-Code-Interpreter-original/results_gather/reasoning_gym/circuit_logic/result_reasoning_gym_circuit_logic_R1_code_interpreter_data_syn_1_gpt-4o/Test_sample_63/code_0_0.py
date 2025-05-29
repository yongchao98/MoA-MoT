# Define the input values
A = 0
B = 0
C = 1
D = 0
E = 1
F = 0

# Define the logic gate functions
def NOT(x):
    return 1 if x == 0 else 0

def NAND(x, y):
    return NOT(x & y)

def AND(x, y):
    return x & y

def XOR(x, y):
    return x ^ y

def NOR(x, y):
    return NOT(x | y)

# Evaluate the circuit step by step
# Negate gates
negate_F = NOT(F)
negate_E = NOT(E)
negate_D = NOT(D)
negate_C = NOT(C)
negate_B = NOT(B)
negate_A = NOT(A)

# NAND gates
nand_1 = NAND(negate_F, negate_E)
nand_2 = NAND(negate_D, negate_C)
nand_3 = NAND(negate_B, negate_A)

# AND gate
and_1 = AND(nand_1, nand_2)

# XOR gates
xor_1 = XOR(and_1, nand_3)
xor_2 = XOR(nand_2, nand_3)
xor_3 = XOR(nand_1, xor_2)

# NOR gate
nor_1 = NOR(xor_1, xor_3)

# Final output
OUT = nor_1

print(OUT)