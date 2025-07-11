def solve_maximal_chromatic_number():
    """
    This script calculates the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    It does so by reasoning about the properties of the graph sum and cycles.
    """

    print("Step-by-step derivation of the maximal chromatic number of G = C_n + C_n + C_n:")
    print("-" * 70)

    # Upper Bound Analysis
    print("\n1. Upper Bound using Product Coloring:")
    print("The chromatic number of a graph sum G = G1 + G2 + ... + Gk is bounded by the product of the individual chromatic numbers.")
    print("chi(G) <= chi(C_n) * chi(C_n) * chi(C_n) = chi(C_n)^3.")
    
    # Chromatic number of a cycle
    chi_cn_even = 2
    chi_cn_odd = 3
    
    upper_bound_even = chi_cn_even ** 3
    upper_bound_odd = chi_cn_odd ** 3
    
    print(f"   - If n is even (n>=4), chi(C_n) = {chi_cn_even}. Thus, chi(G) <= {chi_cn_even}^3 = {upper_bound_even}.")
    print(f"   - If n is odd (n>=3), chi(C_n) = {chi_cn_odd}. Thus, chi(G) <= {chi_cn_odd}^3 = {upper_bound_odd}.")

    # Lower Bound Analysis
    print("\n2. Lower Bound using Independence Number:")
    print("A lower bound for the chromatic number is given by chi(G) >= |V(G)| / alpha(G).")
    print("   - The number of vertices |V(G)| = n^3.")
    print("   - The independence number alpha(G) = alpha(C_n)^3 = (floor(n/2))^3.")
    print("So, chi(G) >= n^3 / (floor(n/2))^3.")
    
    print("\n3. Analyzing Specific Cases:")
    # Case: n is even
    print("   - For even n:")
    print(f"     The lower bound is n^3 / (n/2)^3 = 8.")
    print(f"     Since 8 <= chi(G) <= {upper_bound_even}, the chromatic number is exactly {upper_bound_even}.")

    # Case: n = 3
    n = 3
    print(f"   - For n = {n}:")
    print(f"     C_{n} is the complete graph K_{n}.")
    print(f"     G = K_{n} + K_{n} + K_{n}.")
    print("     The sum of complete graphs is itself a complete graph on the product of the vertex sets.")
    num_vertices = n**3
    print(f"     So, G is the complete graph K_({n}*{n}*{n}) = K_{num_vertices}.")
    chi_G_n3 = num_vertices
    print(f"     The chromatic number of K_{num_vertices} is {chi_G_n3}.")

    # Conclusion
    print("\n4. Conclusion:")
    print("We have found that:")
    print(f"   - For even n, chi(G) = {upper_bound_even}.")
    print(f"   - For odd n, chi(G) <= {upper_bound_odd}.")
    print(f"   - The value of {chi_G_n3} is achieved for n = {n}.")
    
    max_chromatic_number = chi_G_n3
    print("\nTherefore, the maximal chromatic number G can have is determined by the case n=3.")
    print(f"Maximal chromatic number = {max_chromatic_number}")

    print("\nThe final equation for the maximal value (achieved at n=3) is:")
    print(f"chi(K_{{{n}}} + K_{{{n}}} + K_{{{n}}}) = chi(K_{{{n}**3}}) = {n} * {n} * {n} = {n**3}")


if __name__ == '__main__':
    solve_maximal_chromatic_number()
