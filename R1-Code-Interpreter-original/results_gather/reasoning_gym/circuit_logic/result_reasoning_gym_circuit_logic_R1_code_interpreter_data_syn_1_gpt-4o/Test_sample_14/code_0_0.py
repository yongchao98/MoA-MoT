# Given input values
A = 1
B = 0
C = 1
D = 1
E = 1

# Evaluate the circuit step by step
# AND gates
and1 = A and B
and2 = C and D
and3 = E and and1
and4 = and2 and and3

# NOT gate
not1 = not and4

# Final AND gates
and5 = and1 and not1
and6 = and2 and and5
and7 = and3 and and6

# Final output
OUT = and7

print(OUT)