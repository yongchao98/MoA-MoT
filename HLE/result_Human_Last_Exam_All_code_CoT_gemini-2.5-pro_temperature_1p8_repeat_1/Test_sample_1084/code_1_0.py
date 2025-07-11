import sys

def solve():
    """
    Analyzes the six statements about graph properties.
    
    My analysis relies on constructing a counterexample graph family that violates all the given statements.
    Let G be the lexicographic product of a 5-cycle with an independent set of size d/2, denoted C_5[I_{d/2}].
    This graph is d-regular, and has n = 5 * (d/2) = 2.5d vertices.
    
    The number of four-cycle homomorphisms, A, in this graph can be shown to be of the order Theta(d^4).
    The number of C_6' homomorphisms, B, is 0 because the structure of C_5 prevents the necessary connections. A C_6' requires two vertices to be connected by paths of length 1 and 3 simultaneously (in terms of partitions). In C_5, this is not possible.

    For B=0, each statement of the form A = O(B/d^k + T(n,d)) simplifies to A = O(T(n,d)).
    We test each statement with A = Theta(d^4) and n = Theta(d).

    1) A = O(B/d^0.5 + nd^2) => d^4 = O(0 + d*d^2) = O(d^3). False.
    2) A = O(B/d^0.99 + nd^2) => d^4 = O(0 + d*d^2) = O(d^3). False.
    3) A = O(B/d + nd^2) => d^4 = O(0 + d*d^2) = O(d^3). False.
    4) A = O(B/d^1.5 + nd^2) => d^4 = O(0 + d*d^2) = O(d^3). False.
    5) A = O(B/d^1.3 + n d^2.6) => d^4 = O(0 + d*d^2.6) = O(d^3.6). False.
    6) A = O(B/d + nd^1.5) => d^4 = O(0 + d*d^1.5) = O(d^2.5). False.

    Since all statements are false for this family of graphs, they are not true for every graph.
    """
    
    # Based on the step-by-step analysis, none of the statements hold for all d-regular graphs.
    results = ['N', 'N', 'N', 'N', 'N', 'N']
    
    print("".join(results))

solve()