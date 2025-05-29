# Possible values
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Find integer L and M
L, M = None, None
for l in values:
    m = 51 - l
    if m in values and m == 2.4 * l:
        L, M = l, m
        break

# Find integer H and K
H, K = None, None
for k in values:
    h = 1.5 * k
    if h in values and k - h == -8:
        H, K = h, k
        break

# Find integer A and D
A, D = None, None
for a in values:
    d = a - 49
    if d in values:
        A, D = a, d
        break

# Remaining values to assign: E, F, G, I, J, C
remaining_values = [v for v in values if v not in [L, M, H, K, A, D]]

# Find F
F = 3.0 * L
if F not in remaining_values:
    F = None

# Find E and I
E, I = None, None
for e in remaining_values:
    i = 8 - e
    if i in remaining_values:
        E, I = e, i
        break

# Find G
G = 2.0 * I
if G not in remaining_values:
    G = None

# Find C and J
C, J = None, None
for c in remaining_values:
    j = 9 - c
    if j in remaining_values:
        C, J = c, j
        break

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