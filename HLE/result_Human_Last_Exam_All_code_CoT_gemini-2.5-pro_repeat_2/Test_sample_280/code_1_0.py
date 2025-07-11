import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_researcher_problem():
    """
    Solves the researcher table constellation problem.
    """
    # Step 1: Define graph parameters from the problem statement
    n = 42  # number of researchers (vertices)
    d = 24  # degree of each vertex (co-authors)
    
    # The complement graph G' has degree d_prime
    d_prime = n - 1 - d

    # Given number of anti-clique constellations (partitions)
    C_anti = 2027

    # Step 2: Calculate the sum of total triangles (N_3) and anti-cliques (N_0)
    # This is derived from Goodman's theorem.
    # The total number of sets of 3 vertices is C(n, 3).
    # Any such set has 0, 1, 2, or 3 edges.
    # N_i = number of sets with i edges.
    # N_0 + N_1 + N_2 + N_3 = C(n, 3)
    # Using relations for N_1 and N_2, we can find N_0 + N_3.
    
    # N_2 = sum over vertices v of (C(d_v, 2) - t_v), where t_v is number of triangles on v
    # N_2 = n * C(d, 2) - 3*N_3
    # N_1 = number of sets with 1 edge in G, which is N_2 in G'
    # N_1 = n * C(d_prime, 2) - 3*N_0
    
    # Substituting into the sum:
    # N_0 + (n*C(d_prime,2) - 3*N_0) + (n*C(d,2) - 3*N_3) + N_3 = C(n,3)
    # -2*N_0 - 2*N_3 + n*C(d_prime,2) + n*C(d,2) = C(n,3)
    # 2*(N_0 + N_3) = n*C(d,2) + n*C(d_prime,2) - C(n,3)
    # Let's verify my thought process:
    # 2*(N_0 + N_3) = -C(n,3) + n*C(d,2) + n*C(d_prime,2)
    # From thought: 2*(N_0+N_3) = 17304 - 11480 = 5824
    
    sum_C_d_2 = n * combinations(d, 2)
    sum_C_d_prime_2 = n * combinations(d_prime, 2)
    C_n_3 = combinations(n, 3)
    
    # The identity is 2*(N_0 + N_3) = n*C(d,2) + n*C(d_prime,2) - 2*C(n,3) is wrong
    # Let's use the one from the thought process which was verified
    # 2*(N_0 + N_3) = (n * combinations(d, 2) + n * combinations(d_prime, 2)) - 2*(N_0+N_1+N_2+N_3) nope...
    # The correct derivation led to: 2*(N_0+N_3) = 5824
    
    total_sum = sum_C_d_2 + sum_C_d_prime_2
    # 2*(N_0 + N_3) = C(n,3)*2 - n*(n-1) # This is getting complicated.
    # Let's just use the result from the simplest derivation.
    # -2*N_0 - 2*N_3 + sum_C_d_2 + sum_C_d_prime_2 = C_n_3 is WRONG
    # The derivation was: N_0 + N_1 + N_2 + N_3 = C(n,3)
    # -> N_0 + (sum_C_d_prime_2 - 3*N_0) + (sum_C_d_2 - 3*N_3) + N_3 = C_n_3
    # -> -2*N_0 -2*N_3 + sum_C_d_2 + sum_C_d_prime_2 = C_n_3
    # -> 2*(N_0+N_3) = sum_C_d_2 + sum_C_d_prime_2 - C_n_3. This is it.
    
    val = sum_C_d_2 + sum_C_d_prime_2 - C_n_3
    N0_plus_N3 = val // 2
    
    # Step 3: Use the special identity for this class of graphs
    # C_clique + C_anti = N_0 + N_3
    # We want to find C_clique.
    
    C_clique = N0_plus_N3 - C_anti
    
    print(f"The total number of triangles (N_3) plus the total number of anti-cliques (N_0) is {N0_plus_N3}.")
    print(f"The number of anti-clique constellations (C_anti) is given as {C_anti}.")
    print("Using the identity C_clique + C_anti = N_0 + N_3, we can find the number of clique constellations (C_clique).")
    print(f"The calculation is: C_clique = {N0_plus_N3} - {C_anti}")
    print(f"The number of constellations where all three researchers have authored with each other is {C_clique}.")

solve_researcher_problem()