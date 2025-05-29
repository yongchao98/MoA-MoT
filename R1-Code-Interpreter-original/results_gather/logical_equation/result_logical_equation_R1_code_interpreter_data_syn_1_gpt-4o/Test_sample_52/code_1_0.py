# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve for A and H
A = (159 - 141) / 2
H = 159 - A

# Solve for D
D = A + 66

# Solve for I
I = A + 112

# Solve for B, E, and J
# E = 3.0B and E = 2.4J
# Let's try possible values for B and J
possible_B = [num for num in numbers if num * 3 in numbers]
possible_J = [num for num in numbers if num * 2.4 in numbers]

# Find the correct B and J
for b in possible_B:
    e = b * 3
    if e in numbers:
        for j in possible_J:
            if e == j * 2.4:
                B = b
                E = e
                J = j
                break

# Solve for G and F
G = B + 209
F = G - 145

# Assign the remaining numbers to C and K
used_values = {A, B, D, E, F, G, H, I, J}
remaining_values = [num for num in numbers if num not in used_values]

# Assign the remaining values to C and K
C, K = remaining_values

# Create the final result list in alphabetical order
result = [A, B, C, D, E, F, G, H, I, J, K]

print(result)