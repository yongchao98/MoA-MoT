# Define the logic gate functions
def NAND(x, y):
    return int(not (x and y))

def OR(x, y):
    return int(x or y)

def NOT(x):
    return int(not x)

# Input values
A = 1
B = 0
C = 1
D = 0
E = 0
F = 0

# Evaluate the circuit step by step
# First layer of NAND gates
nand1 = NAND(A, B)
nand2 = NAND(C, D)
nand3 = NAND(E, F)

# Second layer of NAND gates
nand4 = NAND(nand1, nand2)
nand5 = NAND(nand3, nand4)

# Third layer of NAND gates
nand6 = NAND(nand5, nand5)

# NOT gates
not1 = NOT(nand6)

# OR gates
or1 = OR(not1, nand5)

# Final output
output = OR(or1, nand4)

# Print the final output
print(output)