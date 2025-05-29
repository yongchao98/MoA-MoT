# Given numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Solve for D and E
E = 6  # From the constraint E = 2.0G, and E - I = -26, E must be 6
D = 9  # From D = 1.5E, D = 1.5 * 6 = 9

# Solve for G
G = 3  # From E = 2.0G, G = E / 2 = 6 / 2 = 3

# Solve for I
I = 32  # From E - I = -26, I = E + 26 = 6 + 26 = 32

# Solve for C
C = 36  # From C + I = 68, C = 68 - I = 68 - 32 = 36

# Solve for A and F
A = 121  # From A + F = 201, and A must be the largest number available
F = 80   # F = 201 - A = 201 - 121 = 80

# Assign remaining numbers to B and H
B = 20  # Remaining numbers are 20 and 120, B is smaller
H = 120 # H is larger than G, so H = 120

# Collect the values in alphabetical order
values = [A, B, C, D, E, F, G, H, I]
print(values)