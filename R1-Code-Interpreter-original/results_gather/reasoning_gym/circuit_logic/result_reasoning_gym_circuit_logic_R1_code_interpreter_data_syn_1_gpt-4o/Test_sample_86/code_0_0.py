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
A, B, C, D, E, F, G = 0, 0, 0, 0, 0, 1, 0

# Evaluate the circuit step by step
# Following the connections from the diagram

# First layer of gates
n1 = NAND(F, G)
n2 = NAND(E, n1)
n3 = NAND(D, n2)
n4 = NAND(C, n3)
n5 = NAND(B, n4)
n6 = NAND(A, n5)

# Second layer of gates
a1 = AND(F, G)
a2 = AND(E, a1)
a3 = AND(D, a2)
a4 = AND(C, a3)
a5 = AND(B, a4)
a6 = AND(A, a5)

# Third layer of gates
n7 = NAND(n6, a6)

# Fourth layer of gates
nor1 = NOR(n7, n5)

# Final output
OUT = NOR(nor1, n4)

# Print the final output
print(OUT)