# Given inputs
A = 1
B = 1

# Negations
not_A = 1 - A
not_B = 1 - B

# AND operations
and1 = A and B
and2 = A and not_B
and3 = not_A and B
and4 = not_A and not_B

# XOR operations
xor1 = and1 ^ and2
xor2 = and3 ^ and4
xor3 = xor1 ^ xor2

# Final output
output = xor3

print(output)