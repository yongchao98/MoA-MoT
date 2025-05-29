# Input values
A, B, C, D, E, F, G, H = 0, 1, 0, 1, 0, 1, 0, 1

# NAND gate function
def NAND(x, y):
    return 0 if x == 1 and y == 1 else 1

# AND gate function
def AND(x, y):
    return 1 if x == 1 and y == 1 else 0

# NOT gate function
def NOT(x):
    return 1 if x == 0 else 0

# OR gate function
def OR(x, y):
    return 1 if x == 1 or y == 1 else 0

# Evaluate the circuit
nand1 = NAND(A, B)
nand2 = NAND(C, D)
nand3 = NAND(E, F)
nand4 = NAND(G, H)

and1 = AND(nand1, nand2)
and2 = AND(nand3, nand4)

not1 = NOT(and1)
not2 = NOT(and2)

output = OR(not1, not2)

print(output)