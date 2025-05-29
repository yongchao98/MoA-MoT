# Given input values
A = 1
B = 0
C = 1
D = 1
E = 0
F = 1
G = 0
H = 0

# Define logic gate functions
def AND(x, y):
    return x & y

def OR(x, y):
    return x | y

def XOR(x, y):
    return x ^ y

def NOT(x):
    return ~x & 1

# Evaluate the circuit step by step
# Start from the innermost gates and work towards the output

# First layer of XOR gates
XOR1 = XOR(H, G)
XOR2 = XOR(XOR1, F)
XOR3 = XOR(XOR2, E)
XOR4 = XOR(XOR3, D)
XOR5 = XOR(XOR4, C)
XOR6 = XOR(XOR5, B)
XOR7 = XOR(XOR6, A)

# Negate the output of some XOR gates
NOT_XOR4 = NOT(XOR4)
NOT_XOR5 = NOT(XOR5)
NOT_XOR6 = NOT(XOR6)

# AND gates
AND1 = AND(NOT_XOR4, NOT_XOR5)
AND2 = AND(AND1, NOT_XOR6)

# OR gates
OR1 = OR(XOR7, AND2)

# Final output
OUT = OR1

# Print the final output
print(OUT)