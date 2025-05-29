# Given input values
A = 1
B = 0
C = 0
D = 0
E = 0
F = 0
G = 0

# Define logic gate functions
def NOT(x):
    return 1 if x == 0 else 0

def AND(x, y):
    return x & y

def XOR(x, y):
    return x ^ y

# Evaluate the circuit step by step
# Start from the innermost gates and move towards the output

# XOR gates
xor1 = XOR(F, G)
xor2 = XOR(E, xor1)
xor3 = XOR(D, xor2)
xor4 = XOR(C, xor3)
xor5 = XOR(B, xor4)
xor6 = XOR(A, xor5)

# NOT gates
not1 = NOT(xor6)
not2 = NOT(xor5)
not3 = NOT(xor4)
not4 = NOT(xor3)
not5 = NOT(xor2)
not6 = NOT(xor1)

# AND gates
and1 = AND(not1, not2)
and2 = AND(not3, not4)
and3 = AND(not5, not6)

# Final AND gate for output
final_output = AND(and1, and2)
final_output = AND(final_output, and3)

# Print the final output
print(final_output)