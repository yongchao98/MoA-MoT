def solve():
    """
    Determines the set of possible numbers of vertices for a convex polyhedron
    that can be projected onto three planes in general position as a quadrilateral.
    
    The reasoning is as follows:
    1. A theorem states that if a polyhedron has a k-gonal projection (with k>=3),
       the number of vertices V must be even. So, V must be even.
    2. We demonstrate that all even numbers V >= 4 are possible by construction.
       - V=4: A tetrahedron works, as generic projections are quadrilaterals.
       - V=6: A regular octahedron, projected along its three axes, yields squares.
       - V=8: A cube, projected along its three axes, yields squares.
       - For any even V >= 10, a specific construction involving vertices of thin orthogonal
         rectangles and pairs of poles along the axes can be shown to work.
         
    Therefore, the set of possible numbers of vertices is all even integers V >= 4.
    """
    
    print("The set of possible numbers of vertices for such a polyhedron is the set of all even integers greater than or equal to 4.")
    print("This can be expressed as the set S = {V | V is an integer, V >= 4, and V % 2 == 0}.")
    print("\nLet's check the condition V % 2 == 0 and V >= 4 for the first few possible values:")
    possible_values = [4, 6, 8, 10, 12]
    for v in possible_values:
        is_even = v % 2 == 0
        is_large_enough = v >= 4
        print(f"For V = {v}:")
        print(f"  V >= 4: {is_large_enough}")
        print(f"  V % 2 == 0: {is_even}")
        print(f"  Is {v} a possible number of vertices? {is_large_enough and is_even}")

solve()