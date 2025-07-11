import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def matchings_in_cycle(n, k):
    """
    Calculates the number of k-matchings in a cycle graph C_n.
    Formula: n/(n-k) * C(n-k, k)
    """
    if k == 0:
        return 1
    if n < 2 * k:
        return 0
    return round(n / (n - k) * combinations(n - k, k))

# Define the parameters for our counterexample
n_vertices = 12
degree = 2
k_matchings = 3

print(f"We will compare two graphs, G1 and G2, that are both bipartite, {degree}-regular, and have {n_vertices} vertices.")
print("-" * 40)

# Graph G1: C12 (A single cycle of 12 vertices)
n1 = 12
g1_matchings = matchings_in_cycle(n1, k_matchings)

print(f"Graph G1 is a cycle graph C{n1}.")
print(f"The number of {k_matchings}-matchings in G1 is calculated using the formula: n/(n-k) * C(n-k, k)")
C_val = combinations(n1 - k_matchings, k_matchings)
print(f"m(C{n1}, {k_matchings}) = {n1}/({n1}-{k_matchings}) * C({n1}-{k_matchings}, {k_matchings})")
print(f"           = {n1}/{n1-k_matchings} * {C_val}")
print(f"           = {g1_matchings}")

print("-" * 40)

# Graph G2: 3C4 (Disjoint union of three C4 cycles)
n_comp = 4
num_comp = 3
print(f"Graph G2 is the disjoint union of {num_comp} cycle graphs C{n_comp}.")
print(f"To find the number of {k_matchings}-matchings in G2, we sum the ways to distribute {k_matchings} edges among the {num_comp} components.")

# First, calculate the number of smaller matchings in a single C4 component
m_c4_0 = matchings_in_cycle(n_comp, 0)
m_c4_1 = matchings_in_cycle(n_comp, 1)
m_c4_2 = matchings_in_cycle(n_comp, 2)
m_c4_3 = matchings_in_cycle(n_comp, 3)

print(f"\nNumber of matchings in a single C4:")
print(f"m(C4, 0) = {m_c4_0}")
print(f"m(C4, 1) = {m_c4_1} (the number of edges)")
print(f"m(C4, 2) = {m_c4_2}")
print(f"m(C4, 3) = {m_c4_3} (a C4 is too small for a 3-matching)")

# A 3-matching in G2 can be formed by partitioning the 3 edges among the three C4s.
# The possible partitions of 3 are (1,1,1) and (2,1,0).

# Case 1: Partition (1,1,1) -> 1 edge from each C4.
ways_111 = m_c4_1 * m_c4_1 * m_c4_1
print("\nCase 1: Partition (1,1,1). One edge is chosen from each of the three C4 components.")
print(f"Number of ways = m(C4,1) * m(C4,1) * m(C4,1) = {m_c4_1} * {m_c4_1} * {m_c4_1} = {ways_111}")

# Case 2: Partition (2,1,0) -> 2 edges from one C4, 1 from another, 0 from the third.
# There are 3 choices for the component that contributes 2 edges.
ways_210 = combinations(num_comp, 1) * (m_c4_2 * m_c4_1 * m_c4_0)
print("\nCase 2: Partition (2,1,0). Two edges from one C4, one from another, and zero from the third.")
print(f"Number of ways = C(3,1) * (m(C4,2) * m(C4,1) * m(C4,0)) = 3 * ({m_c4_2} * {m_c4_1} * {m_c4_0}) = {ways_210}")

g2_matchings = ways_111 + ways_210
print(f"\nThe total number of {k_matchings}-matchings in G2 is the sum of these cases:")
print(f"Total = {ways_111} + {ways_210} = {g2_matchings}")

print("-" * 40)
print(f"Number of 3-matchings in G1 (C12) = {g1_matchings}")
print(f"Number of 3-matchings in G2 (3C4) = {g2_matchings}")

if g1_matchings != g2_matchings:
    print("\nSince the numbers are different, two bipartite, d-regular graphs on n vertices do not necessarily have the same number of 3-matchings.")
else:
    print("\nFor this example, the numbers are the same. This is not a counterexample.")
