# Input assignments
A = 1
B = 1
C = 1
D = 0
E = 1
F = 1
G = 0
H = 0
I = 0

# Define logic gate functions
def NAND(x, y):
    return int(not (x and y))

def XOR(x, y):
    return int(x != y)

def OR(x, y):
    return int(x or y)

def NOT(x):
    return int(not x)

# Evaluate the circuit step by step
# Following the circuit from inputs to output

# First layer of NAND gates
nand1 = NAND(A, B)
nand2 = NAND(C, D)
nand3 = NAND(E, F)
nand4 = NAND(G, H)
nand5 = NAND(I, nand1)

# Second layer of NAND gates
nand6 = NAND(nand2, nand3)
nand7 = NAND(nand4, nand5)

# XOR gates
xor1 = XOR(nand6, nand7)
xor2 = XOR(nand3, nand4)

# OR gates
or1 = OR(xor1, xor2)
or2 = OR(nand5, nand6)

# Final output
OUT = OR(or1, or2)

# Print the final output
print(OUT)