# Input assignments
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
    return all(args)

def NAND(*args):
    return not all(args)

def NOR(*args):
    return not any(args)

def XOR(a, b):
    return a != b

def NOT(a):
    return not a

# Evaluate the circuit step by step
# Layer 1
nand1 = NAND(H, G)
nand2 = NAND(F, E)
nand3 = NAND(D, C)
nand4 = NAND(B, A)

# Layer 2
nand5 = NAND(nand1, nand2)
nand6 = NAND(nand3, nand4)

# Layer 3
nor1 = NOR(nand5, nand6)

# Layer 4
nand7 = NAND(nor1, nand5)

# Final output
OUT = NOR(nand7, nand6)

# Print the final output
print(int(OUT))