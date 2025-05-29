from itertools import permutations

numbers = [5, 41, 38, 81, 14]

# Try all permutations of the numbers
for perm in permutations(numbers):
    # Try different combinations of operations
    # (a * b) + (c * d) + e
    if perm[0] * perm[1] + perm[2] * perm[3] + perm[4] == 450:
        print(f"{perm[0]} * {perm[1]} + {perm[2]} * {perm[3]} + {perm[4]}")
    # (a * b) + (c * d) - e
    if perm[0] * perm[1] + perm[2] * perm[3] - perm[4] == 450:
        print(f"{perm[0]} * {perm[1]} + {perm[2]} * {perm[3]} - {perm[4]}")
    # (a * b) - (c * d) + e
    if perm[0] * perm[1] - perm[2] * perm[3] + perm[4] == 450:
        print(f"{perm[0]} * {perm[1]} - {perm[2]} * {perm[3]} + {perm[4]}")
    # (a * b) - (c * d) - e
    if perm[0] * perm[1] - perm[2] * perm[3] - perm[4] == 450:
        print(f"{perm[0]} * {perm[1]} - {perm[2]} * {perm[3]} - {perm[4]}")
    # (a + b) * (c + d) + e
    if (perm[0] + perm[1]) * (perm[2] + perm[3]) + perm[4] == 450:
        print(f"({perm[0]} + {perm[1]}) * ({perm[2]} + {perm[3]}) + {perm[4]}")
    # (a + b) * (c + d) - e
    if (perm[0] + perm[1]) * (perm[2] + perm[3]) - perm[4] == 450:
        print(f"({perm[0]} + {perm[1]}) * ({perm[2]} + {perm[3]}) - {perm[4]}")