# Assign the numbers to the letters based on the derived relationships
A = 121
B = 3
C = 4
D = 80
E = 9
F = 1
G = 36
H = 16
I = 27

# Verify the inequalities
assert C - E == -5
assert D - F == 79
assert G == 4 * E
assert H == 4 * C
assert E - G == -27
assert E + H == 25
assert A + H == 137

# Output the results in alphabetical order
result = [A, B, C, D, E, F, G, H, I]
print(result)