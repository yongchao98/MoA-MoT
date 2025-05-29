# Define logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def XOR(x, y):
    return x ^ y

def NOR(x, y):
    return ~(x | y) & 1

def Negate(x):
    return ~x & 1

# Inputs
A = 0
B = 0
C = 1

# Evaluate the circuit
and1 = AND(A, B)
and2 = AND(B, C)

nand1 = NAND(A, B)
nand2 = NAND(B, C)
nand3 = NAND(A, C)
nand4 = NAND(and1, and2)

xor1 = XOR(A, B)
xor2 = XOR(nand1, nand2)

negate = Negate(nand3)

nor1 = NOR(nand4, negate)
nor2 = NOR(xor2, nor1)

# Output
print(nor2)