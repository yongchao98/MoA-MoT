# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Function to check if a set of values satisfies the conditions
def check_conditions(A, B, C, D, E, F, G, H, I):
    return (
        E + G == 127 and
        H == 4.0 * D and
        A + D == 140 and
        B == 1.6 * D and
        A - G == -1 and
        G > B and
        C + D == 29 and
        G - C == 112
    )

# Iterate over all permutations of the values
from itertools import permutations

for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    if check_conditions(A, B, C, D, E, F, G, H, I):
        print([A, B, C, D, E, F, G, H, I])
        break