def demonstrate_contradiction():
    """
    This function explains the logical contradiction arising from the problem's conditions.
    It doesn't search for a graph, but rather proves that no such graph can exist.
    """
    print("Analyzing the properties of the graph G with n vertices:")
    print("-" * 50)

    # From Property 3: The graph contains exactly n copies of C5 (cycles of length 5).
    # From Property 4: No three of these C5s can share a common vertex.

    print("Let's use a double-counting argument on the (vertex, C5) pairs.")
    print("A (vertex, C5) pair means the vertex is part of that C5 cycle.")

    # 1. Counting by summing over all C5 cycles:
    print("\nStep 1: Counting by summing over the C5 cycles.")
    print("The graph has 'n' copies of C5.")
    print("Each C5 has 5 vertices.")
    print("Total number of (vertex, C5) incidences = (Number of C5s) * (Vertices per C5)")
    print("Total incidences = n * 5")
    total_incidences_by_cycles = "5 * n"

    # 2. Counting by summing over all vertices:
    print("\nStep 2: Counting by summing over the vertices.")
    print("Let c(v) be the number of C5s that vertex 'v' belongs to.")
    print("Property 4 states that no three C5s share a vertex, so for any vertex 'v', c(v) <= 2.")
    print("This means a vertex can be in 0, 1, or 2 cycles.")
    
    print("\nLet n_k be the number of vertices belonging to exactly k cycles.")
    print("The total number of vertices is n, so we have our first equation:")
    print("n_0 + n_1 + n_2 = n  (Equation 1)")

    print("\nThe total number of (vertex, C5) incidences is the sum of c(v) over all vertices:")
    print("Total incidences = (0 * n_0) + (1 * n_1) + (2 * n_2)")
    total_incidences_by_vertices = "n_1 + 2 * n_2"
    
    # 3. Equating the two counts:
    print("\nStep 3: Equating the two counts to form our second equation.")
    print(f"From Step 1 and 2: {total_incidences_by_vertices} = {total_incidences_by_cycles}")
    print("n_1 + 2 * n_2 = 5 * n  (Equation 2)")

    # 4. Solving the system of equations:
    print("\nStep 4: Analyzing the system of equations.")
    print("We have:")
    print("1) n_0 + n_1 + n_2 = n")
    print("2) n_1 + 2 * n_2 = 5 * n")

    print("\nLet's express n_1 from Equation 1: n_1 = n - n_0 - n_2")
    print("Substitute this into Equation 2:")
    print("(n - n_0 - n_2) + 2 * n_2 = 5 * n")
    print("n - n_0 + n_2 = 5 * n")
    print("n_2 - n_0 = 4 * n")

    print("\nNow, let's use n_2 = n_0 + 4*n and substitute it back into Equation 1:")
    print("n_0 + n_1 + (n_0 + 4 * n) = n")
    print("2 * n_0 + n_1 + 4 * n = n")
    print("This simplifies to the final equation:")
    final_equation_lhs_n0 = 2
    final_equation_lhs_n1 = 1
    final_equation_rhs_n = -3
    print(f"{final_equation_lhs_n0} * n_0 + {final_equation_lhs_n1} * n_1 = {final_equation_rhs_n} * n")
    
    # 5. The contradiction:
    print("\nStep 5: The Contradiction.")
    print("The variables n_0 and n_1 represent counts of vertices, so they must be non-negative integers (>= 0).")
    print("The number of vertices 'n' for a composite number must be positive (> 0).")
    print(f"Therefore, the left side ({final_equation_lhs_n0}*n_0 + {final_equation_lhs_n1}*n_1) must be a non-negative number.")
    print(f"The right side ({final_equation_rhs_n}*n) must be a negative number.")
    print("A non-negative number cannot be equal to a negative number.")

    print("-" * 50)
    print("Conclusion: The given properties are self-contradictory. No such graph exists for any positive n.")

if __name__ == '__main__':
    demonstrate_contradiction()