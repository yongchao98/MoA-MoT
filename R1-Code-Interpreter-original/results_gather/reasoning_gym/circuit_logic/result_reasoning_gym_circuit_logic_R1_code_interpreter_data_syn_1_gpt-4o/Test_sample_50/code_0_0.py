# Define logic gate functions
def NOT(x):
    return 1 if x == 0 else 0

def AND(x, y):
    return x & y

def NAND(x, y):
    return NOT(AND(x, y))

def XOR(x, y):
    return x ^ y

# Inputs
A = 1
B = 1
C = 1
D = 1

# Negate gates
not_A = NOT(A)
not_B = NOT(B)
not_C = NOT(C)
not_D = NOT(D)

# XOR gates
xor1 = XOR(C, D)
xor2 = XOR(B, xor1)
xor3 = XOR(A, xor2)

# AND gates
and1 = AND(not_C, xor3)
and2 = AND(not_B, xor2)
and3 = AND(not_A, and2)
and4 = AND(and1, and3)

# NAND gates
nand1 = NAND(and4, not_D)
nand2 = NAND(and3, nand1)

# Final AND gate
OUT = AND(and4, nand2)

# Output the result
print(OUT)