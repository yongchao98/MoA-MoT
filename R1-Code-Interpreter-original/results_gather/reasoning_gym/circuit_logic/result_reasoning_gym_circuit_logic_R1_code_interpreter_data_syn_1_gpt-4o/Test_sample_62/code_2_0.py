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
and1 = AND(H, G)
and2 = AND(F, E)
and3 = AND(D, C)
and4 = AND(B, A)

# Layer 2
nand1 = NAND(and1, and2)
nand2 = NAND(and3, and4)

# Layer 3
nor1 = NOR(nand1, nand2)

# Layer 4
nand3 = NAND(nor1, nand1)

# Layer 5
nor2 = NOR(nand3, nand2)

# Final Output
OUT = nor2

print(OUT)