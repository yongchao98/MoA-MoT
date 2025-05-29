# Input values
A = 0
B = 0
C = 0
D = 1
E = 1
F = 1
G = 1

# Define logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def XOR(x, y):
    return x ^ y

def NOT(x):
    return ~x & 1

# Evaluate the circuit
# Layer 1
nG = NOT(G)
nF = NOT(F)
nE = NOT(E)
nD = NOT(D)
nC = NOT(C)
nB = NOT(B)
nA = NOT(A)

# Layer 2
and1 = AND(nG, nF)
and2 = AND(nE, nD)
and3 = AND(nC, nB)
and4 = AND(nA, G)

# Layer 3
nand1 = NAND(and1, and2)
nand2 = NAND(and3, and4)

# Layer 4
xor1 = XOR(nand1, nand2)

# Layer 5
nand3 = NAND(xor1, F)

# Layer 6
xor2 = XOR(nand3, E)

# Final Output
OUT = xor2

print(OUT)