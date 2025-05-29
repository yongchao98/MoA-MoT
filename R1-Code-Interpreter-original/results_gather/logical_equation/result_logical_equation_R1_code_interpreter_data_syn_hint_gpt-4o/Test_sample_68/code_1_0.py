# Assign known values based on reasoning
B = 28
L = 7
J = 15
M = 45
C = 24
G = 3
K = 16
E = 5
I = 36
H = 10

# Remaining values to assign: A, D, F
remaining_values = [1, 2, 50]

# Check the constraint J > D
# J = 15, so D must be less than 15
# Possible values for D are 1, 2

# Assign D = 2 (since it's less than 15)
D = 2

# Now assign the remaining values to A and F
A, F = 1, 50

# Print the final assignment
result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
print(result)