# Define the logic gate functions
def AND(x, y):
    return x & y

def NAND(x, y):
    return ~(x & y) & 1

def NOR(x, y):
    return ~(x | y) & 1

def NOT(x):
    return ~x & 1

# Input values
A = 1
B = 1
C = 0
D = 1

# Calculate the outputs of the gates
nand1 = NAND(A, B)
nand2 = NAND(C, D)
nand3 = NAND(A, C)
nand4 = NAND(B, D)

and1 = AND(nand1, nand2)
and2 = AND(nand3, nand4)

negate = NOT(and1)

final_output = NOR(negate, and2)

# Print the final output
print(final_output)