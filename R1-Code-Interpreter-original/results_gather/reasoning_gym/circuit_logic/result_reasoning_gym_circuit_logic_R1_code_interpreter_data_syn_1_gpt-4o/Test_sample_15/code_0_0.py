# Define the logic gate functions
def nand_gate(x, y):
    return int(not (x and y))

def negate(x):
    return int(not x)

def or_gate(x, y):
    return int(x or y)

# Inputs
A = 0
B = 1
C = 0
D = 1
E = 1
F = 1
G = 0
H = 1

# Evaluate the circuit step by step
# Following the circuit diagram from top to bottom, left to right

# First layer of NAND gates
nand1 = nand_gate(H, G)
nand2 = nand_gate(nand1, F)
nand3 = nand_gate(nand2, E)
nand4 = nand_gate(nand3, D)
nand5 = nand_gate(nand4, C)
nand6 = nand_gate(nand5, B)
nand7 = nand_gate(nand6, A)

# Negate the output of the last NAND gate
negate1 = negate(nand7)

# OR gate with the negated output
output = or_gate(negate1, nand7)

# Print the final output
print(output)