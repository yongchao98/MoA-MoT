# Define the logic gate functions
def AND(*args):
    return all(args)

def NAND(*args):
    return not all(args)

def XOR(a, b):
    return a != b

def NOR(*args):
    return not any(args)

def NOT(a):
    return not a

# Given inputs
A = 1
B = 0
C = 1
D = 0
E = 0
F = 1
G = 1
H = 1

# Evaluate the circuit step by step
# Layer 1
and1 = AND(H, G)
and2 = AND(F, E)
and3 = AND(D, C)
and4 = AND(B, A)

# Layer 2
and5 = AND(and1, and2)
and6 = AND(and3, and4)

# Layer 3
nor1 = NOR(and5, and6)

# Layer 4
nand1 = NAND(and5, and6)

# Layer 5
xor1 = XOR(nor1, nand1)

# Final Output
output = NOR(xor1, nand1)

# Print the final output
print(output)