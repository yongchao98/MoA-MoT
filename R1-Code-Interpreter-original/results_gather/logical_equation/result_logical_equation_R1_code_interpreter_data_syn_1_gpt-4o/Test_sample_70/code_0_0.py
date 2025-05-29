# Possible values
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Solve for L and M
L = 51 / 3.4
M = 2.4 * L

# Solve for H and K
K = 16
H = 1.5 * K

# Solve for A and D
A = 50
D = 1

# Remaining values to assign: E, F, G, I, J, C
remaining_values = [v for v in values if v not in [L, M, H, K, A, D]]

# Solve for F
F = 3.0 * L

# Solve for E and I
E = 5
I = 3

# Solve for G
G = 2.0 * I

# Solve for C and J
C = 5
J = 4

# Assign values to letters
letter_values = {
    'A': A,
    'B': None,
    'C': C,
    'D': D,
    'E': E,
    'F': F,
    'G': G,
    'H': H,
    'I': I,
    'J': J,
    'K': K,
    'L': L,
    'M': M
}

# Find the remaining value for B
remaining_values = [v for v in values if v not in letter_values.values()]
letter_values['B'] = remaining_values[0]

# Create the final list in alphabetical order
final_values = [letter_values[letter] for letter in sorted(letter_values.keys())]

print(final_values)