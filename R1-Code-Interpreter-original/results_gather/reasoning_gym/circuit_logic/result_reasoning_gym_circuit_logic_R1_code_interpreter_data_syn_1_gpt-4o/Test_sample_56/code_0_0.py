# Define the logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def NOR(x, y):
    return ~(x | y) & 1

def NOT(x):
    return ~x & 1

# Given inputs
A = 1
B = 0
C = 0
D = 0
E = 1

# Evaluate the circuit step by step
# First layer of NAND gates
nand1 = NAND(A, B)
nand2 = NAND(C, D)
nand3 = NAND(E, nand1)

# Second layer of NAND gates
nand4 = NAND(nand2, nand3)

# AND gates
and1 = AND(nand1, nand2)
and2 = AND(nand3, nand4)

# Final NOR gate
output = NOR(and1, and2)

# Print the final output
print(output)