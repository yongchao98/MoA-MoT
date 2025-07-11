import math

def combinations(n, k):
    """
    Helper function to calculate the number of combinations (n choose k).
    """
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# --- Introduction ---
print("This script calculates the symmetry factors for the two second-order vacuum bubble diagrams in phi^4 theory and sums them.")
print("The symmetry factor S is calculated using the formula: S = (V! * (n!)^V) / C")
print("Here, V=2 (vertices), n=4 (from phi^4), and C is the number of ways to form the diagram's contractions.\n")

# --- Diagram 1: The "Setting-Sun" Diagram ---
print("1. The 'Setting-Sun' Diagram (4 propagators between 2 vertices)")

# Number of contractions (C_sun):
# This is the number of ways to connect the 4 legs of vertex 1 to the 4 legs of vertex 2.
# This is equivalent to the number of permutations of the 4 legs on one side, which is 4!.
c_sun = math.factorial(4)

# Symmetry factor (S_sun):
V = 2
n = 4
s_sun = (math.factorial(V) * (math.factorial(n))**V) / c_sun

print(f"   - Number of ways to form (C_sun): 4! = {c_sun}")
print(f"   - Symmetry factor (S_sun): (2! * (4!)^2) / {c_sun} = {int(s_sun)}\n")


# --- Diagram 2: The "Figure-Eight" Diagram ---
print("2. The 'Figure-Eight' Diagram (1 loop on each vertex, 2 propagators between them)")

# Number of contractions (C_fig8):
# Ways to form a loop on vertex 1: C(4, 2)
# Ways to form a loop on vertex 2: C(4, 2)
# Ways to connect the remaining 2 legs from vertex 1 to vertex 2: 2
c42 = combinations(4, 2)
c_fig8 = c42 * c42 * 2

# Symmetry factor (S_fig8):
s_fig8 = (math.factorial(V) * (math.factorial(n))**V) / c_fig8

print(f"   - Number of ways to form (C_fig8): C(4,2) * C(4,2) * 2 = {c42} * {c42} * 2 = {c_fig8}")
print(f"   - Symmetry factor (S_fig8): (2! * (4!)^2) / {c_fig8} = {int(s_fig8)}\n")


# --- Final Summation ---
print("--- Total Sum ---")
total_s = s_sun + s_fig8
# The final equation is printed below as requested.
print(f"{int(s_sun)} + {int(s_fig8)} = {int(total_s)}")
<<<64>>>