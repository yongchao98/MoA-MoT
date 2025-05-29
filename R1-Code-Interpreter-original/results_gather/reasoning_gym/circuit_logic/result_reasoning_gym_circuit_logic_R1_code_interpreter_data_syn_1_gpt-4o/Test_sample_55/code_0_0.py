# Define the logic gate functions
def NOT(x):
    return 1 if x == 0 else 0

def AND(x, y):
    return x & y

def NAND(x, y):
    return NOT(AND(x, y))

def XOR(x, y):
    return x ^ y

def NOR(x, y):
    return NOT(x | y)

# Input values
A = 1
B = 1
C = 1
D = 1

# Negate gates
neg_A = NOT(A)
neg_B = NOT(B)
neg_C = NOT(C)
neg_D = NOT(D)

# XOR gates
xor1 = XOR(neg_A, neg_B)
xor2 = XOR(neg_C, neg_D)
xor3 = XOR(xor1, xor2)

# AND gates
and1 = AND(neg_A, neg_B)
and2 = AND(neg_C, neg_D)
and3 = AND(and1, and2)

# NAND gates
nand1 = NAND(xor1, xor2)
nand2 = NAND(and1, and2)

# NOR gates
nor1 = NOR(xor3, nand1)
nor2 = NOR(and3, nand2)

# Final output
output = nor2

print(output)