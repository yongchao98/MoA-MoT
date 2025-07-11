def analyze_graph_problem():
    """
    Analyzes the properties of the graph described in the problem
    and demonstrates the logical contradiction that arises from them.
    """
    print("This script analyzes the properties of a graph G with n vertices to find the smallest composite n.")
    
    print("\nLet's break down the given properties:")
    print("1. G is 7-regular: Every vertex has a degree of 7.")
    print("2. Chromatic number χ(G) = 5: G cannot be colored with 4 colors.")
    print("3. G contains exactly n copies of C5 (cycles of length 5).")
    print("4. No three of these C5s can share a common vertex.")
    
    print("\nNow, let's turn these properties into mathematical statements.")
    
    print("\nLet S be the set of all 5-cycles in G. Property 3 means |S| = n.")
    print("Let c(v) be the number of cycles in S that pass through a vertex v.")
    
    print("\nFrom Property 4, 'No three of these C5s can share a common vertex', we can deduce a constraint on c(v).")
    print("If c(v) were 3 or more, it would mean that vertex v is part of at least three cycles from S.")
    print("This would violate the condition. Therefore, c(v) must be less than 3 for every vertex v.")
    print("This gives us the inequality: c(v) <= 2 for all v in G.")

    print("\nWe will now use a double-counting argument on the (vertex, cycle) pairs.")
    print("We count the total number of times a vertex appears in a cycle in two ways.")

    print("\nMethod 1: Summing over all cycles in S.")
    print("Each cycle is a C5, meaning it has 5 vertices. There are n cycles in S.")
    print("Total incidences = (Number of cycles) * (Vertices per cycle)")
    print("Total incidences = n * 5")
    
    print("\nMethod 2: Summing over all vertices in G.")
    print("The total number of incidences is the sum of c(v) over all n vertices.")
    print("Total incidences = sum(c(v) for all v in G)")
    
    print("\nFrom our analysis of Property 4, we know that c(v) <= 2 for every vertex.")
    print("So, we can find an upper bound for this sum:")
    print("sum(c(v)) <= (Number of vertices) * (Max c(v))")
    print("sum(c(v)) <= n * 2")
    
    print("\nNow, we combine the findings.")
    print("From Method 1, the total number of incidences is exactly 5*n.")
    print("From Method 2, the total number of incidences is at most 2*n.")
    
    print("\nThis leads to the following inequality:")
    # The final equation as requested by the user.
    print("5 * n <= 2 * n")
    
    print("\nThis inequality simplifies to 3 * n <= 0, which means n <= 0.")
    print("However, the number of vertices 'n' in a graph must be a positive integer.")
    print("This is a logical contradiction, which proves that no graph can satisfy all the given properties for any positive n.")
    
    print("\nConclusion: The problem statement is contradictory. No such graph exists, so there is no smallest composite n for which it exists.")

if __name__ == '__main__':
    analyze_graph_problem()