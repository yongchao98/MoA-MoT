# This script determines the matching between airfoil shapes and their pressure plots.

# The matching is based on the following analysis:
# 1. Airfoils A, B, C, D are thin. Airfoils E, F, G, H are thick.
# 2. Pressure plots 1, 2, 4, 6 have sharp suction peaks characteristic of thin airfoils.
# 3. Pressure plots 3, 5, 7, 8 have rounded peaks characteristic of thick airfoils.
# 4. Within each group, lift (area of the Cp loop) increases with camber.

# Ranking by camber/lift:
# Thin: A > B > C > D(symmetric)
# Thick: E > F > G > H(symmetric)

# Ranking by lift/area from plots:
# Thin-airfoil plots: Area(6) > Area(1) > Area(2) > Area(4)
# Thick-airfoil plots: Area(5) > Area(7) > Area(3) > Area(8)

# Assembling the pairs (Airfoil -> Plot):
A = 6  # Thin, highest camber -> Highest lift in thin group
B = 1  # Thin, high camber -> High lift in thin group
C = 2  # Thin, low camber -> Low lift in thin group
D = 4  # Thin, symmetric -> Lowest lift in thin group
E = 5  # Thick, highest camber -> Highest lift in thick group
F = 7  # Thick, high camber -> High lift in thick group
G = 3  # Thick, low camber -> Low lift in thick group
H = 8  # Thick, symmetric -> Lowest lift in thick group

# The final answer is the sequence of plot numbers for airfoils A through H.
final_sequence = f"{A}{B}{C}{D}{E}{F}{G}{H}"

print(final_sequence)