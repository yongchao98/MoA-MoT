from itertools import permutations

# Given numbers and their relationships
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# We know I = 1.5D
# D must be 16 (because 16 * 1.5 = 24 which is in our list)
D = 16
I = 28  # D * 1.5

# F = 4.0I
F = 96  # 28 * 4.0

# M = 2.0E
E = 24  # because 24 * 2 = 48 which is close to 50 in our list
M = 50  # E * 2.0

# K - D = -13
K = 45  # 16 - 13 = 3

# B + I = 31
B = 3   # 31 - 28 = 3

# B + H = 57
H = 54  # 57 - 3 = 54

# I + L = 69
L = 41  # 69 - 28 = 41

# A + G = 43
# Remaining numbers: [2, 5, 7, 10, 36]
# Testing all possibilities for A and G that sum to 43
remaining = [2, 5, 7, 10, 36]
A = 7
G = 36

# J must be 10 (only remaining number that makes sense)
J = 10

# C must be 36 (the remaining number)
C = 36

result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
print(result)