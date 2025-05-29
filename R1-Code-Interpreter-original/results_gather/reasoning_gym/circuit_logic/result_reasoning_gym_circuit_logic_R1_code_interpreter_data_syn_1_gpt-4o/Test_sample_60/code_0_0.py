# Define the logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def XOR(x, y):
    return x ^ y

def NOT(x):
    return ~x & 1

# Input values
A = 0
B = 1
C = 0
D = 1
E = 1

# Evaluate the circuit step by step
# First layer of gates
xor1 = XOR(A, B)
xor2 = XOR(C, D)
xor3 = XOR(E, xor1)

# Second layer of gates
not1 = NOT(xor2)
xor4 = XOR(not1, xor3)

# Third layer of gates
and1 = AND(xor2, xor3)
nand1 = NAND(xor1, xor4)

# Fourth layer of gates
xor5 = XOR(and1, nand1)

# Final output
OUT = xor5

# Print the final output
print(OUT)