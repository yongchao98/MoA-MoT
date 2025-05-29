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
and3 = and1 and E
and4 = and2 and (not E)
and5 = and3 and and4
and6 = and5 and and2
and7 = and6 and and1

# Final output
OUT = and7

print(OUT)