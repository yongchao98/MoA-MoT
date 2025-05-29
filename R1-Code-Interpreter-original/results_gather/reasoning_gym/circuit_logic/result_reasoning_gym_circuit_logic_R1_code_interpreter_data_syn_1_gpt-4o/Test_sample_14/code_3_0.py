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
and3 = and2 and E
and4 = and1 and (not and3)
and5 = and4 and and2
and6 = and5 and and1

# Final output
OUT = and6

print(OUT)