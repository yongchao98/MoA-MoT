# Possible values for each letter
values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Check each combination to find a consistent solution
from itertools import permutations

for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    if (B + E == 170 and
        H - C == -46 and
        K == 3.0 * G and
        B + I == 230 and
        C - H == 46 and
        D == 3.0 * H and
        C - D == -104 and
        H - J == 36 and
        A == 2.4 * E and
        I - E == 60 and
        C - K == 112):
        result = [A, B, C, D, E, F, G, H, I, J, K]
        break

print(result)