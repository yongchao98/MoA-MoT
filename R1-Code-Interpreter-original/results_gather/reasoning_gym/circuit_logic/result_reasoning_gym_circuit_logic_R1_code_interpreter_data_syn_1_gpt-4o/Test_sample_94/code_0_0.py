# Define the logic gate functions
def XOR(a, b):
    return a ^ b

def NAND(a, b):
    return not (a and b)

def Negate(a):
    return not a

def OR(a, b):
    return a or b

# Input values
A, B, C, D, E, F, G, H, I = 1, 1, 1, 1, 1, 1, 0, 1, 0

# Evaluate the circuit
XOR1 = XOR(A, B)
XOR2 = XOR(C, D)
XOR3 = XOR(E, F)
XOR4 = XOR(G, H)
XOR5 = XOR(I, XOR1)
XOR6 = XOR(XOR2, XOR3)
XOR7 = XOR(XOR4, XOR5)
XOR8 = XOR(XOR6, XOR7)

NAND1 = NAND(XOR8, XOR5)

Negate1 = Negate(XOR8)
Negate2 = Negate(XOR7)
Negate3 = Negate(XOR6)

OR1 = OR(Negate1, NAND1)
OR2 = OR(Negate2, Negate3)
OR3 = OR(OR1, OR2)

# Final output
print(int(OR3))