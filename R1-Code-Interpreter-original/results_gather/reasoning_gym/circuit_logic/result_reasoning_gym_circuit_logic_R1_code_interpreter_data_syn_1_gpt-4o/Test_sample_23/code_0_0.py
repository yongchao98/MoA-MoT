# Define logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def XOR(x, y):
    return x ^ y

def NOR(x, y):
    return ~(x | y) & 1

def NOT(x):
    return ~x & 1

# Inputs
A = 0
B = 0

# Evaluate the circuit
and1 = AND(A, B)
nand1 = NAND(A, B)
xor1 = XOR(A, B)

negate1 = NOT(and1)
negate2 = NOT(xor1)

nor1 = NOR(nand1, negate1)
nor2 = NOR(nand1, negate2)

# Final output
output = nor2

print(output)