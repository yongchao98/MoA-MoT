# Define logic gate functions
def AND(*inputs):
    return all(inputs)

def NAND(*inputs):
    return not all(inputs)

def NOR(*inputs):
    return not any(inputs)

def NOT(input):
    return not input

# Given input values
A = 1
B = 1
C = 0
D = 1
E = 0
F = 1
G = 0
H = 0
I = 1

# Evaluate the circuit step by step
# NAND gates
nand_1 = NAND(I, H)
nand_2 = NAND(nand_1, G)
nand_3 = NAND(nand_2, F)
nand_4 = NAND(nand_3, E)
nand_5 = NAND(nand_4, D)
nand_6 = NAND(nand_5, C)
nand_7 = NAND(nand_6, B)
nand_8 = NAND(nand_7, A)

# AND gates
and_1 = AND(nand_8, nand_7)
and_2 = AND(and_1, nand_6)
and_3 = AND(and_2, nand_5)
and_4 = AND(and_3, nand_4)
and_5 = AND(and_4, nand_3)
and_6 = AND(and_5, nand_2)
and_7 = AND(and_6, nand_1)

# NOR gate for final output
output = NOR(and_7, nand_8)

# Print the final output
print(int(output))