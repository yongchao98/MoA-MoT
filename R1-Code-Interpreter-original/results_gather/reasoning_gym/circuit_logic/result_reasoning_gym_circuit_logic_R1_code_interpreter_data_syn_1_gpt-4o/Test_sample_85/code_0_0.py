# Define logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def XOR(x, y):
    return x ^ y

def NOT(x):
    return ~x & 1

def NOR(x, y):
    return ~(x | y) & 1

# Given inputs
A = 0
B = 0
C = 0

# Evaluate the circuit step by step
# First layer of gates
and1 = AND(A, B)
and2 = AND(B, C)
and3 = AND(A, C)

# Second layer of gates
xor1 = XOR(and1, and2)
xor2 = XOR(and2, and3)
xor3 = XOR(and1, and3)

# Third layer of gates
not_xor1 = NOT(xor1)
not_xor2 = NOT(xor2)
not_xor3 = NOT(xor3)

# Fourth layer of gates
nor1 = NOR(not_xor1, not_xor2)
nor2 = NOR(not_xor2, not_xor3)
nor3 = NOR(not_xor1, not_xor3)

# Final output
output = NOR(nor1, nor2)

print(output)