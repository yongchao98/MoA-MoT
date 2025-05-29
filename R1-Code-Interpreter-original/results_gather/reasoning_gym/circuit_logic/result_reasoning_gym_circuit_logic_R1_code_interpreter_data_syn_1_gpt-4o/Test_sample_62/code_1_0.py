# Input values
A = 0
B = 1
C = 0
D = 0
E = 0
F = 1
G = 0
H = 0

# Define logic gate functions
def AND(*args):
    return int(all(args))

def NAND(*args):
    return int(not all(args))

def NOR(*args):
    return int(not any(args))

def NOT(x):
    return int(not x)

# Evaluate the circuit step by step
# Layer 1
nand1 = NAND(H, G)
nand2 = NAND(F, E)
nand3 = NAND(D, C)
nand4 = NAND(B, A)

# Layer 2
nor1 = NOR(nand1, nand2)
nor2 = NOR(nand3, nand4)

# Layer 3
nand5 = NAND(nor1, nor2)

# Layer 4
nor3 = NOR(nand5, nor1)

# Final Output
OUT = nor3

print(OUT)