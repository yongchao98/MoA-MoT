# Define the logic gates
def xor_gate(*inputs):
    return sum(inputs) % 2

def and_gate(*inputs):
    return int(all(inputs))

def not_gate(input):
    return int(not input)

# Given inputs
A = 0
B = 1
C = 0
D = 0
E = 1

# Evaluate the circuit
# Start from the innermost gates and work towards the output

# First layer of XOR gates
xor1 = xor_gate(A, B)
xor2 = xor_gate(C, D)
xor3 = xor_gate(E, xor1)

# NOT gate
not1 = not_gate(xor2)

# Second layer of XOR gates
xor4 = xor_gate(not1, xor3)

# AND gates
and1 = and_gate(xor1, xor2)
and2 = and_gate(xor3, xor4)

# Final XOR gate for output
output = xor_gate(and1, and2)

# Print the final output
print(output)