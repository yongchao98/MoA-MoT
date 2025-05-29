# Define the input values
A = 1
B = 0
C = 1
D = 0
E = 0
F = 1
G = 0

# Define the logic gate functions
def NOT(x):
    return 1 if x == 0 else 0

def AND(*args):
    return 1 if all(args) else 0

def NAND(*args):
    return NOT(AND(*args))

def XOR(*args):
    return 1 if sum(args) % 2 == 1 else 0

# Evaluate the circuit step by step
# Negate gates
nA = NOT(A)
nB = NOT(B)
nC = NOT(C)
nD = NOT(D)
nE = NOT(E)
nF = NOT(F)
nG = NOT(G)

# XOR gates
xor1 = XOR(nF, nG)
xor2 = XOR(nE, xor1)
xor3 = XOR(nD, xor2)
xor4 = XOR(nC, xor3)
xor5 = XOR(nB, xor4)
xor6 = XOR(nA, xor5)

# NAND gates
nand1 = NAND(nF, nG)
nand2 = NAND(nE, nand1)
nand3 = NAND(nD, nand2)
nand4 = NAND(nC, nand3)
nand5 = NAND(nB, nand4)
nand6 = NAND(nA, nand5)

# AND gates
and1 = AND(xor6, nand6)

# Final output
OUT = and1

# Print the final output
print(OUT)