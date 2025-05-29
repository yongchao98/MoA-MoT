# Define the logic gate functions
def AND(*args):
    return all(args)

def NAND(*args):
    return not all(args)

def XOR(a, b):
    return a != b

def NOT(a):
    return not a

# Given input values
A, B, C, D, E, F, G, H, I, J, K = 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0

# Evaluate the circuit step by step
# Layer 1: Negations
not_F = NOT(F)
not_I = NOT(I)
not_J = NOT(J)
not_K = NOT(K)

# Layer 2: AND gates
and_1 = AND(not_K, J)
and_2 = AND(not_J, I)
and_3 = AND(not_I, H)
and_4 = AND(not_F, G)
and_5 = AND(E, D)
and_6 = AND(C, B)
and_7 = AND(A, not_K)

# Layer 3: NAND gates
nand_1 = NAND(and_1, and_2)
nand_2 = NAND(and_3, and_4)
nand_3 = NAND(and_5, and_6)

# Layer 4: XOR gates
xor_1 = XOR(nand_1, nand_2)
xor_2 = XOR(nand_3, and_7)

# Final XOR gate for output
output = XOR(xor_1, xor_2)

# Print the final output
print(output)