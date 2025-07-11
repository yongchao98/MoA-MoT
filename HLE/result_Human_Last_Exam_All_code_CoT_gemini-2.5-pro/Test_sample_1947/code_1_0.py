# This script calculates and prints the coefficients for the number of
# closed tree-like walks of length 6 in a simple graph X.

# The formula is of the form:
# c_1*e + c_2*k + c_3*p + c_4*sum(deg(v) choose 2) + c_5*sum(deg(v) choose 3)

# Based on combinatorial counting, we determined the following coefficients:

# c_1 corresponds to walks on a P_2 (single edge) traversed 6 times.
c_1 = 2

# c_2 corresponds to K_3 subgraphs. Tree-like walks cannot have a cycle as
# their underlying structure, so this coefficient is 0.
c_2 = 0

# c_3 corresponds to walks on a P_4 where each of the 3 edges is traversed twice.
c_3 = 6

# c_4 corresponds to walks on a P_3 where one edge is traversed twice and the other four times.
c_4 = 14

# c_5 corresponds to walks on a K_{1,3} where each of the 3 edges is traversed twice.
c_5 = 12

print(f"c_1 = {c_1}")
print(f"c_2 = {c_2}")
print(f"c_3 = {c_3}")
print(f"c_4 = {c_4}")
print(f"c_5 = {c_5}")