# Define the logic gate functions
def AND(*args):
    return all(args)

def NAND(*args):
    return not all(args)

def XOR(a, b):
    return a != b

def NOT(a):
    return not a

# Input values
A, B, C, D, E, F, G, H, I, J = 1, 0, 1, 1, 1, 1, 0, 0, 1, 1

# Evaluate the circuit
# Start from the innermost gates and work outwards

# First layer of gates
and1 = AND(J, I)
and2 = AND(and1, H)
and3 = AND(and2, G)
and4 = AND(and3, F)
and5 = AND(and4, E)
and6 = AND(and5, D)
and7 = AND(and6, C)
and8 = AND(and7, B)
and9 = AND(and8, A)

# Second layer of gates
nand1 = NAND(and9, and8)
nand2 = NAND(and7, and6)

# Third layer of gates
xor1 = XOR(nand1, nand2)
xor2 = XOR(and5, and4)

# Fourth layer of gates
not1 = NOT(xor2)

# Fifth layer of gates
xor3 = XOR(not1, xor1)

# Final output
output = xor3

# Print the final output
print(output)