# Define the logic gate functions
def NAND(x, y):
    return int(not (x and y))

def AND(x, y):
    return x and y

def NOR(x, y):
    return int(not (x or y))

def NOT(x):
    return int(not x)

# Given inputs
A, B, C, D, E, F, G, H, I = 1, 1, 0, 1, 0, 1, 0, 0, 1

# Evaluate the circuit step by step
# Layer 1
nand1 = NAND(I, H)
nand2 = NAND(nand1, G)
nand3 = NAND(nand2, F)
nand4 = NAND(nand3, E)
nand5 = NAND(nand4, D)
nand6 = NAND(nand5, C)
nand7 = NAND(nand6, B)
nand8 = NAND(nand7, A)

# Layer 2
and1 = AND(nand8, nand7)
and2 = AND(and1, nand6)
and3 = AND(and2, nand5)
and4 = AND(and3, nand4)
and5 = AND(and4, nand3)
and6 = AND(and5, nand2)
and7 = AND(and6, nand1)

# Layer 3
nor1 = NOR(and7, and6)
nor2 = NOR(nor1, and5)
nor3 = NOR(nor2, and4)
nor4 = NOR(nor3, and3)
nor5 = NOR(nor4, and2)
nor6 = NOR(nor5, and1)

# Final output
OUT = nor6

# Print the final output
print(OUT)