# Input values
A, B, C, D, E, F = 1, 1, 1, 1, 1, 0

# Negated values
A_, B_, C_, D_, E_, F_ = 0, 0, 0, 0, 0, 1

# Gate evaluations
and1 = A and B and C and D and E and F
and2 = A_ and B_ and C_ and D_ and E_ and F_
and3 = A and B and C and D and E
and4 = A_ and B_ and C_ and D_ and E_
nand1 = not (A and B and C and D and E)
nand2 = not (A_ and B_ and C_ and D_ and E_)

# Final output
output = and1 and and2 and and3 and and4 and nand1 and nand2
print(output)