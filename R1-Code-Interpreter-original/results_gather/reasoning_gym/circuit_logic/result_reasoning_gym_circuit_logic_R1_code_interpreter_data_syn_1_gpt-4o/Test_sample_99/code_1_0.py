# Define the logic gate functions
def xor_gate(a, b):
    return a ^ b

def nand_gate(a, b):
    return not (a and b)

def not_gate(a):
    return not a

# Given input values
A = 1
B = 1
C = 1
D = 0
E = 0
F = 0
G = 0
H = 0

# Evaluate the circuit step by step
# From the circuit diagram, we need to evaluate the gates in order

# First layer of gates
xor1 = xor_gate(A, B)
xor2 = xor_gate(C, D)
xor3 = xor_gate(E, F)
xor4 = xor_gate(G, H)

# Second layer of gates
xor5 = xor_gate(xor1, xor2)
xor6 = xor_gate(xor3, xor4)

# Third layer of gates
nand1 = nand_gate(xor5, xor6)

# Fourth layer of gates
not1 = not_gate(nand1)

# Final XOR gate for output
final_output = xor_gate(not1, xor6)

# Print the final output
print(int(final_output))