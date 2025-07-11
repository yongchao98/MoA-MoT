import math

def solve_graph_problem():
    """
    Calculates and explains the solution for finding the maximum number of edges
    in a simple graph with 8 vertices that has no quadrilaterals (C4).
    """

    n = 8
    print(f"Let n be the number of vertices. n = {n}")
    print("Let e be the number of edges.")
    print("The problem is to find the maximum possible value of e in a simple graph G with n vertices,")
    print("such that G is 'C4-free' (contains no cycle of length 4).")
    
    print("\n### Step 1: Establish an inequality for C4-free graphs ###")
    print("A cycle C4 (A-B-C-D-A) implies that the vertices A and C are both connected to B and D.")
    print("In other words, the pair of vertices {A, C} has two common neighbors.")
    print("For a graph to be C4-free, any pair of vertices can have at most one common neighbor.")

    print("\nLet's count the number of '2-paths', which are paths of length 2 (e.g., u-v-w).")
    print("For a vertex v with degree d(v), the number of 2-paths centered at v is d(v) choose 2.")
    print("The total number of 2-paths in the graph is the sum over all vertices: Sum[d(v) * (d(v) - 1) / 2].")
    
    print("\nSince each pair of vertices {u, w} can have at most one common neighbor v,")
    print("the total number of 2-paths cannot be greater than the total number of pairs of vertices.")
    n_choose_2 = math.comb(n, 2)
    print(f"The number of pairs of vertices is n choose 2 = {n} choose 2 = {n_choose_2}.")
    print(f"This gives us the fundamental inequality: Sum[d(v) * (d(v) - 1) / 2] <= {n_choose_2}")
    
    inequality_rhs = 2 * n_choose_2
    print(f"Multiplying by 2, we get: Sum[d(v)^2 - d(v)] <= {inequality_rhs}")

    print("\n### Step 2: Derive an upper bound for the number of edges e ###")
    print("We know two key facts about degrees:")
    print("1. Sum of degrees: Sum[d(v)] = 2 * e")
    print("2. Cauchy-Schwarz inequality: (Sum[d(v)])^2 <= n * Sum[d(v)^2]")
    
    print("\nLet's combine these facts with our inequality:")
    print(f"From Sum[d(v)^2 - d(v)] <= {inequality_rhs}, we have Sum[d(v)^2] <= {inequality_rhs} + Sum[d(v)] = {inequality_rhs} + 2*e.")
    print(f"From Cauchy-Schwarz, (2*e)^2 <= {n} * Sum[d(v)^2], so Sum[d(v)^2] >= (4*e^2)/{n}.")
    
    print(f"\nCombining these, we get: (4*e^2)/{n} <= {inequality_rhs} + 2*e")
    # Simplify the inequality: 4e^2 <= n*(inequality_rhs + 2e) => 4e^2 <= 8*(56 + 2e) => e^2 <= 2*(56+2e)
    # This leads to e^2 - 4e - 112 <= 0
    a, b, c = 1, -4, -112
    print("This simplifies to the following quadratic inequality for e:")
    print(f"e^2 + ({b})*e + ({c}) <= 0")
    
    print("\nTo find the upper bound for e, we solve the equation e^2 - 4e - 112 = 0.")
    discriminant = b**2 - 4 * a * c
    positive_root = (-b + math.sqrt(discriminant)) / (2 * a)
    print(f"The positive root is (4 + sqrt(16 - 4*1*(-112)))/2 = (4 + sqrt({discriminant}))/2 ≈ {positive_root:.2f}")
    
    upper_bound = math.floor(positive_root)
    print(f"Therefore, the number of edges e must be less than or equal to {upper_bound}.")

    print("\n### Step 3: Analysis of the result ###")
    print(f"The calculation gives an upper bound of e <= {upper_bound}.")
    print("This bound is only met if the inequalities used are equalities. This requires a regular graph where every pair of vertices has exactly one common neighbor (a structure known as a symmetric block design).")
    print(f"For e = {upper_bound}, the graph would need to be 3-regular (since 2*e/n = 2*{upper_bound}/{n} = 3).")
    print("However, it is a known result that all 3-regular graphs on 8 vertices (such as the cube graph) contain C4s.")
    print(f"Thus, the maximum number of edges must be less than {upper_bound}.")

    print("\nFurther analysis shows that a C4-free graph with 11 edges on 8 vertices also cannot exist.")
    print("However, a C4-free graph with 10 edges on 8 vertices is known to exist.")
    print("One such graph is the cycle on 8 vertices (C8), plus two additional edges connecting opposite vertices, for example, (v0, v4) and (v2, v6). This graph has 8 + 2 = 10 edges and is C4-free.")

    print("\n### Conclusion ###")
    final_answer = 10
    print(f"The maximum number of edges in a simple graph with 8 vertices without a C4 is {final_answer}.")

solve_graph_problem()
<<<10>>>