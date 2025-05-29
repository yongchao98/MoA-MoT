# Define the logic gates as functions
def AND(*inputs):
    return int(all(inputs))

def NAND(*inputs):
    return int(not all(inputs))

def XOR(a, b):
    return int(a != b)

def OR(*inputs):
    return int(any(inputs))

def NOT(a):
    return int(not a)

# Given inputs
A = 0
B = 0
C = 0
D = 1
E = 0

# Evaluate the circuit step by step
# Negate gates
not_A = NOT(A)
not_B = NOT(B)
not_C = NOT(C)
not_D = NOT(D)
not_E = NOT(E)

# First layer of gates
and1 = AND(A, B)
and2 = AND(C, D)
nand1 = NAND(and1, and2)

# Second layer of gates
and3 = AND(nand1, E)
xor1 = XOR(and3, not_E)

# Third layer of gates
and4 = AND(not_C, xor1)
and5 = AND(not_B, and4)

# Fourth layer of gates
nand2 = NAND(and5, not_A)
or1 = OR(nand2, not_D)

# Final output
output = OR(or1, not_E)

# Print the final output
print(output)