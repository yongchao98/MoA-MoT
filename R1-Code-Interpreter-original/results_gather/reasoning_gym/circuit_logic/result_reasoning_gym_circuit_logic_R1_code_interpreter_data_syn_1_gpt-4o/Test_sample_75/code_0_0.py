# Define logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def NOT(x):
    return ~x & 1

# Given inputs
A = 1
B = 1
C = 1
D = 1

# Evaluate the circuit
# First layer of gates
nand1 = NAND(A, B)
nand2 = NAND(C, D)
nand3 = NAND(B, C)
nand4 = NAND(A, D)

# Second layer of gates
not1 = NOT(nand1)
not2 = NOT(nand2)
not3 = NOT(nand3)
not4 = NOT(nand4)

# Third layer of gates
and1 = AND(not1, not2)
and2 = AND(not3, not4)

# Final output
OUT = AND(and1, and2)

print(OUT)