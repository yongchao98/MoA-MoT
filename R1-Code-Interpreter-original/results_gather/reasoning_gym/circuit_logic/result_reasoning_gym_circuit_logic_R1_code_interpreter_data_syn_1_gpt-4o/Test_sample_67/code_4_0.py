# Input values
A = 0
B = 1
C = 0
D = 1

# Negate gates
not_A = 1 - A
not_B = 1 - B
not_C = 1 - C
not_D = 1 - D

# Evaluate the circuit
# First layer of AND gates
and1 = A and B
and2 = C and D
and3 = not_A and not_B
and4 = not_C and not_D

# XOR gates
xor1 = and1 ^ and2
xor2 = and3 ^ and4

# Second layer of AND gates
and5 = xor1 and xor2

# Final AND gate with negated inputs
final_and = and5 and (not_A and not_B) and (not_C and not_D)

# Output
print(final_and)