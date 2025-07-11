import numpy as np

def count_inversions(p):
    """Counts the number of inversions in a permutation."""
    n = len(p)
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inversions += 1
    return inversions

def invert_permutation(p):
    """Computes the inverse of a permutation."""
    # Add 1 to make it 1-based index, easier to invert
    p_1based = [x + 1 for x in p]
    inverse_p_1based = [0] * len(p_1based)
    for i, val in enumerate(p_1based):
        # Invert: if p[i] = val, then inv_p[val] = i
        inverse_p_1based[val - 1] = i + 1
    # Subtract 1 to return to 0-based index
    return [x - 1 for x in inverse_p_1based]

# --- Main Calculation ---
# Grid size
n = 5

# The permutation sigma from the X coordinates (1-based: [4, 5, 1, 2, 3])
# We use 0-based indexing for arrays in Python: [3, 4, 0, 1, 2]
sigma = [3, 4, 0, 1, 2]

# 1. Calculate Writhe
# First, find the inverse of sigma
sigma_inv = invert_permutation(sigma)

# Count inversions for sigma and its inverse
inv_sigma = count_inversions(sigma)
inv_sigma_inv = count_inversions(sigma_inv)

# Writhe is the sum of inversions
writhe = inv_sigma + inv_sigma_inv

# 2. Calculate Number of Right Cusps
# All n 'X' markers are right cusps
right_cusps_X = n

# Count right cusps at 'O' markers
right_cusps_O = 0
for i in range(n):
    # i is 0-indexed, but represents column i+1
    # sigma[i] is 0-indexed, but represents row sigma(i+1)
    vertical_dir = (sigma[i] - i)
    horizontal_dir = ((i + 1) - (sigma_inv[i] + 1))
    if vertical_dir * horizontal_dir < 0:
        right_cusps_O += 1
        
total_right_cusps = right_cusps_X + right_cusps_O

# 3. Calculate Thurston-Bennequin number
tb = writhe - total_right_cusps

# Print the results step-by-step
print(f"The permutation sigma is: {[x + 1 for x in sigma]}")
print(f"The inverse permutation sigma_inv is: {[x + 1 for x in sigma_inv]}")
print(f"Number of inversions in sigma: {inv_sigma}")
print(f"Number of inversions in sigma_inv: {inv_sigma_inv}")
print(f"The writhe of the diagram is the number of crossings = {inv_sigma} + {inv_sigma_inv} = {writhe}")
print(f"Number of right cusps from 'X' markings: {right_cusps_X}")
print(f"Number of right cusps from 'O' markings: {right_cusps_O}")
print(f"Total number of right cusps = {right_cusps_X} + {right_cusps_O} = {total_right_cusps}")
print("\nThe maximal Thurston-Bennequin number is calculated as:")
print(f"tb = writhe - total_right_cusps")
print(f"tb = {writhe} - {total_right_cusps} = {tb}")

print(f"\nFinal Answer: {tb}")
<<<3>>>