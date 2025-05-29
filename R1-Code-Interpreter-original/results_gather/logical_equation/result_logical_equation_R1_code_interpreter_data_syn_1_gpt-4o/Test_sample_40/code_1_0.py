# Remaining values after assigning K, D, C, E, I, J
remaining_values = [9, 16, 36, 48, 75, 80]

# Find possible values for A, F, B, and H
for A in remaining_values:
    F = 52 - A
    B = 3 * F
    H = A + 44
    if F in remaining_values and B in remaining_values and H in remaining_values:
        print(f"A: {A}, F: {F}, B: {B}, H: {H}")