# Define the logic gate functions
def XOR(a, b):
    return a ^ b

def AND(a, b):
    return a & b

def NAND(a, b):
    return ~(a & b) & 1

def NOR(a, b):
    return ~(a | b) & 1

# Inputs
A = 0
B = 1
C = 0
D = 1
E = 1
F = 0

# Evaluate the circuit
XOR1 = XOR(E, F)
XOR2 = XOR(D, XOR1)
XOR3 = XOR(C, XOR2)
XOR4 = XOR(B, XOR3)
XOR5 = XOR(A, XOR4)

AND1 = AND(C, D)
AND2 = AND(AND1, E)
AND3 = AND(AND2, F)

NAND1 = NAND(B, C)
NAND2 = NAND(NAND1, D)
NAND3 = NAND(NAND2, E)

NOR1 = NOR(AND3, NAND3)
NOR2 = NOR(XOR5, NOR1)

# Final output
print(NOR2)