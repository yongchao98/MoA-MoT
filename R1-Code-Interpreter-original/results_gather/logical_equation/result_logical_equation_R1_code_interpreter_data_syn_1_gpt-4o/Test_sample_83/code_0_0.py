# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solving I - J = 51 and I + J = 141
I = (141 + 51) // 2
J = (141 - 51) // 2

# Solving H + I = 103 and H - I = -89
H = (103 - 89) // 2
I = (103 + 89) // 2

# Solving E = 4.0H
E = 4 * H

# Solving G - E = -13
G = E - 13

# Solving G = 3.0A
A = G // 3

# Solving E = 2.8M
M = E // 2.8

# Assigning values to letters
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
values = [A, None, None, None, E, None, G, H, I, J, None, None, M]

# Finding remaining values
remaining_numbers = [num for num in numbers if num not in values]

# Assign remaining numbers to B, C, D, F, K, L
values[1] = remaining_numbers[0]  # B
values[2] = remaining_numbers[1]  # C
values[3] = remaining_numbers[2]  # D
values[5] = remaining_numbers[3]  # F
values[10] = remaining_numbers[4] # K
values[11] = remaining_numbers[5] # L

print(values)