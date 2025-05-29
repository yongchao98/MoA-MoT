# Define the available numbers
available_numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve for G and K using G + K = 375 and G = 1.5K
K = 375 / 2.5
G = 1.5 * K

# Solve for F using F - G = -177
F = G - 177

# Solve for D using F = 2.4D
D = F / 2.4

# Solve for B using F = 3.0B
B = F / 3.0

# Solve for A using A + D = 29
A = 29 - D

# Solve for I using I - D = 19
I = D + 19

# Solve for H using A = 3.0H
H = A / 3.0

# Solve for J using J - F = 73
J = F + 73

# Determine C using G > C
C = min([num for num in available_numbers if num < G])

# Assign the values to the letters
values = [A, B, C, D, E, F, G, H, I, J, K]

# Print the values
print(values)