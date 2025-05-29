# Define logic gate functions
def XOR(a, b):
    return a ^ b

def AND(a, b):
    return a & b

def NAND(a, b):
    return ~(a & b) & 1

def NOR(a, b):
    return ~(a | b) & 1

def NOT(a):
    return ~a & 1

# Inputs
A = 1
B = 1
C = 0

# Evaluate gates
XOR1 = XOR(A, B)
XOR2 = XOR(B, C)
AND1 = AND(A, B)
AND2 = AND(B, C)
NAND1 = NAND(A, B)
NAND2 = NAND(B, C)

NOR1 = NOR(XOR1, XOR2)
NOR2 = NOR(AND1, AND2)
NOR3 = NOR(NAND1, NAND2)

Negate_XOR2 = NOT(XOR2)

# Final output
OUT = NOR(NOR1, NOR2)
OUT = NOR(OUT, NOR3)
OUT = NOR(OUT, Negate_XOR2)

print(OUT)