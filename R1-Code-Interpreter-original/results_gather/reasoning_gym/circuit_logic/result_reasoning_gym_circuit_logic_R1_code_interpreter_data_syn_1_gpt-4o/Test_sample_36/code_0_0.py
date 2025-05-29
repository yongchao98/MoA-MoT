# Input assignments
A, B, C, D, E, F, G, H, I, J, K = 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1

# Define logic gate functions
def XOR(a, b):
    return a ^ b

def NAND(a, b):
    return not (a and b)

def AND(a, b):
    return a and b

def NOT(a):
    return not a

# Evaluate the circuit step by step
# Layer 1
x1 = XOR(J, K)
x2 = XOR(I, x1)
x3 = XOR(H, x2)
x4 = XOR(G, x3)
x5 = XOR(F, x4)
x6 = XOR(E, x5)
x7 = XOR(D, x6)
x8 = XOR(C, x7)
x9 = XOR(B, x8)
x10 = XOR(A, x9)

# Layer 2
n1 = NOT(x10)
n2 = NOT(x9)
n3 = NOT(x8)
n4 = NOT(x7)
n5 = NOT(x6)
n6 = NOT(x5)
n7 = NOT(x4)
n8 = NOT(x3)
n9 = NOT(x2)
n10 = NOT(x1)

# Layer 3
nand1 = NAND(n1, n2)
nand2 = NAND(n3, n4)
nand3 = NAND(n5, n6)
nand4 = NAND(n7, n8)
nand5 = NAND(n9, n10)

# Layer 4
and1 = AND(nand1, nand2)
and2 = AND(nand3, nand4)
and3 = AND(and1, and2)
and4 = AND(and3, nand5)

# Final output
OUT = and4

print(int(OUT))