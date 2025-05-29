# Input values
A = 0
B = 1
C = 1
D = 0
E = 0
F = 1
G = 1
H = 0
I = 1
J = 1

# Negation function
def negate(x):
    return 1 if x == 0 else 0

# XOR function
def xor(x, y):
    return x ^ y

# AND function
def and_gate(*args):
    return int(all(args))

# NAND function
def nand(x, y):
    return negate(and_gate(x, y))

# Evaluate the circuit
# Start from the bottom of the circuit and work upwards

# First layer of XOR gates
xor1 = xor(I, J)
xor2 = xor(H, xor1)
xor3 = xor(G, xor2)
xor4 = xor(F, xor3)
xor5 = xor(E, xor4)
xor6 = xor(D, xor5)
xor7 = xor(C, xor6)
xor8 = xor(B, xor7)
xor9 = xor(A, xor8)

# NAND gates
nand1 = nand(xor8, xor9)
nand2 = nand(xor7, xor8)
nand3 = nand(xor6, xor7)
nand4 = nand(xor5, xor6)
nand5 = nand(xor4, xor5)
nand6 = nand(xor3, xor4)
nand7 = nand(xor2, xor3)
nand8 = nand(xor1, xor2)

# Negation gates
neg1 = negate(nand8)
neg2 = negate(nand7)
neg3 = negate(nand6)
neg4 = negate(nand5)
neg5 = negate(nand4)
neg6 = negate(nand3)
neg7 = negate(nand2)
neg8 = negate(nand1)

# AND gates
and1 = and_gate(neg1, neg2)
and2 = and_gate(neg3, neg4)
and3 = and_gate(neg5, neg6)
and4 = and_gate(neg7, neg8)

# Final AND gate for output
output = and_gate(and1, and2, and3, and4)

# Print the final output
print(output)