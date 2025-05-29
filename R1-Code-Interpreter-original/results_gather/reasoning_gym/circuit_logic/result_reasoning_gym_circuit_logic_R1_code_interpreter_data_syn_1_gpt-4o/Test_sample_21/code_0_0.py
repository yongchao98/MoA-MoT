# Define the logic gate functions
def AND(*inputs):
    return int(all(inputs))

def NAND(*inputs):
    return int(not all(inputs))

def NOR(*inputs):
    return int(not any(inputs))

def NOT(input):
    return int(not input)

# Given inputs
A, B, C, D, E, F, G, H, I = 0, 0, 0, 0, 1, 1, 1, 1, 0

# Evaluate the circuit step by step
# Following the circuit diagram from inputs to output

# First layer of gates
nand1 = NAND(I, H)
nand2 = NAND(nand1, G)
nand3 = NAND(nand2, F)
nand4 = NAND(nand3, E)
nand5 = NAND(nand4, D)
nand6 = NAND(nand5, C)
nand7 = NAND(nand6, B)
nand8 = NAND(nand7, A)

# Second layer of gates
and1 = AND(nand8, nand7)
and2 = AND(nand6, nand5)
and3 = AND(nand4, nand3)
and4 = AND(nand2, nand1)

# Third layer of gates
nor1 = NOR(and1, and2)
nor2 = NOR(and3, and4)

# Final output
OUT = NOR(nor1, nor2)

# Print the final output
print(OUT)