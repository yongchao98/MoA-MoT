# Define the logic gate functions
def AND(*args):
    return all(args)

def NAND(*args):
    return not all(args)

def XOR(a, b):
    return a != b

def NOT(a):
    return not a

def NOR(*args):
    return not any(args)

# Given inputs
A = 0
B = 1
C = 0
D = 0

# Evaluate the circuit step by step
# First layer of AND gates
and1 = AND(A, B)
and2 = AND(C, D)

# XOR gates
xor1 = XOR(and1, C)
xor2 = XOR(and2, B)

# NOR gate
nor1 = NOR(xor1, xor2)

# NOT gate
not1 = NOT(D)

# Second layer of AND gates
and3 = AND(nor1, not1)

# Final NOR gate
output = NOR(and3, A)

# Print the final output
print(int(output))