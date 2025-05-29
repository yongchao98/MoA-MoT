# Define the input values
A = 0
B = 1
C = 0
D = 1
E = 0
F = 1
G = 1
H = 0
I = 0
J = 0

# Define the logic gate functions
def AND(x, y):
    return x & y

def OR(x, y):
    return x | y

def XOR(x, y):
    return x ^ y

def NAND(x, y):
    return ~(x & y) & 1

def NOT(x):
    return ~x & 1

# Evaluate the circuit
# Start from the bottom of the circuit and work upwards

# XOR gates
xor1 = XOR(I, J)
xor2 = XOR(H, xor1)
xor3 = XOR(G, xor2)
xor4 = XOR(F, xor3)
xor5 = XOR(E, xor4)
xor6 = XOR(D, xor5)
xor7 = XOR(C, xor6)
xor8 = XOR(B, xor7)
xor9 = XOR(A, xor8)

# NAND gates
nand1 = NAND(xor5, xor6)
nand2 = NAND(xor4, nand1)
nand3 = NAND(xor3, nand2)

# AND gates
and1 = AND(nand2, nand3)
and2 = AND(xor2, and1)

# OR gates
or1 = OR(and2, xor9)
or2 = OR(nand3, or1)

# Final output
output = or2

# Print the final output
print(output)