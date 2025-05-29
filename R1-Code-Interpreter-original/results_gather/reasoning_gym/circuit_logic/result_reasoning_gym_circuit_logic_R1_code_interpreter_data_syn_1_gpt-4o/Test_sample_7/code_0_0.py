# Define the logic gate functions
def NOT(x):
    return 1 if x == 0 else 0

def NAND(x, y):
    return 0 if x == 1 and y == 1 else 1

def XOR(x, y):
    return 1 if x != y else 0

# Given inputs
A = 1
B = 1
C = 1
D = 1

# Evaluate the circuit step by step
# First layer of NAND gates
nand1 = NAND(A, B)
nand2 = NAND(C, D)

# Second layer of NAND gates
nand3 = NAND(nand1, nand2)
nand4 = NAND(nand1, C)
nand5 = NAND(nand2, D)

# NOT gate
not1 = NOT(nand3)

# XOR gates
xor1 = XOR(nand4, nand5)
xor2 = XOR(not1, xor1)

# Final output
OUT = xor2

print(OUT)