# Given inputs
A = 1
B = 0
C = 1
D = 0
E = 1
F = 1
G = 0

# Negation function
def negate(x):
    return 1 if x == 0 else 0

# XOR function
def xor(x, y):
    return x ^ y

# NOR function
def nor(x, y):
    return 1 if (x == 0 and y == 0) else 0

# Evaluate the circuit
# Step 1: Evaluate XOR gates
xor1 = xor(F, G)  # XOR between F and G
xor2 = xor(E, xor1)  # XOR between E and result of xor1
xor3 = xor(D, xor2)  # XOR between D and result of xor2
xor4 = xor(C, xor3)  # XOR between C and result of xor3
xor5 = xor(B, xor4)  # XOR between B and result of xor4
xor6 = xor(A, xor5)  # XOR between A and result of xor5

# Step 2: Evaluate negation gates
neg1 = negate(xor3)  # Negate the result of xor3
neg2 = negate(xor5)  # Negate the result of xor5

# Step 3: Evaluate NOR gates
nor1 = nor(neg1, xor6)  # NOR between neg1 and xor6
nor2 = nor(neg2, nor1)  # NOR between neg2 and nor1

# Final output
OUT = nor2

# Print the final output
print(OUT)