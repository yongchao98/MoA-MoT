# Define the logic gates
def nand_gate(*inputs):
    return 0 if all(inputs) else 1

def xor_gate(*inputs):
    return sum(inputs) % 2

def not_gate(input):
    return 1 - input

# Inputs
A = 0
B = 1
C = 1
D = 0
E = 1

# Evaluate the circuit step by step
# Following the circuit diagram from left to right, top to bottom

# First layer of NAND gates
nand1 = nand_gate(A, B)
nand2 = nand_gate(C, D)
nand3 = nand_gate(E, nand1)

# Second layer of NAND gates
nand4 = nand_gate(nand2, nand3)

# Negate outputs
neg1 = not_gate(nand1)
neg2 = not_gate(nand2)
neg3 = not_gate(nand3)
neg4 = not_gate(nand4)

# XOR gates
xor1 = xor_gate(neg1, neg2)
xor2 = xor_gate(neg3, neg4)

# Final XOR gate for output
output = xor_gate(xor1, xor2)

# Print the final output
print(output)