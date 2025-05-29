# Define the logic gate functions
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

# Input values
A = 1
B = 1
C = 0
D = 0
E = 0
F = 1
G = 0

# Evaluate the circuit step by step
# Start from the innermost gates and work towards the output

# First layer of gates
nand1 = NAND(F, G)
nand2 = NAND(E, nand1)
nand3 = NAND(D, nand2)
nand4 = NAND(C, nand3)
nand5 = NAND(B, nand4)
nand6 = NAND(A, nand5)

xor1 = XOR(F, G)
xor2 = XOR(E, xor1)
xor3 = XOR(D, xor2)
xor4 = XOR(C, xor3)
xor5 = XOR(B, xor4)
xor6 = XOR(A, xor5)

# Second layer of gates
nor1 = NOR(nand6, xor6)

# Third layer of gates
and1 = AND(nor1, nand6)

# Final output
output = and1

# Print the final output
print(output)