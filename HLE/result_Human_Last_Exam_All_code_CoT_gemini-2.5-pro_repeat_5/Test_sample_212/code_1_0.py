import numpy as np

def demonstrate_actions_on_torus():
    """
    This function demonstrates that for the torus T^2, the two actions of its
    fundamental group on the fiber of its universal cover are the same.
    """
    
    # The fundamental group of the torus, pi_1(T^2), is isomorphic to Z^2.
    # The fiber over a point in the universal cover is also Z^2.
    # We can represent elements of the group and points in the fiber as 2D integer vectors.
    
    # Let g be an element of the fundamental group pi_1(T^2).
    # We choose an example g = (m, n), where m=2, n=3.
    g = np.array([2, 3])
    
    # Let h be a point in the fiber p^{-1}(x_0), also identified with an element of Z^2.
    # We choose an example h = (k, l), where k=5, l=7.
    h = np.array([5, 7])
    
    print("Let g be an element of the fundamental group pi_1(T^2) = Z^2.")
    print(f"g = ({g[0]}, {g[1]})")
    print("\nLet h be a point in the fiber, also identified with an element of Z^2.")
    print(f"h = ({h[0]}, {h[1]})")
    print("-" * 40)

    # Action 1: Holonomy around loops.
    # This action corresponds to right multiplication in the group (h * g).
    # For the additive group Z^2, this is calculated as h + g.
    result_holonomy = h + g
    
    print("Action 1: Holonomy (corresponds to Right Multiplication)")
    print("The action of g on h is calculated as h + g:")
    print(f"   ({h[0]}, {h[1]}) + ({g[0]}, {g[1]}) = ({result_holonomy[0]}, {result_holonomy[1]})")
    print("-" * 40)

    # Action 2: Deck transformations.
    # This action corresponds to left multiplication in the group (g * h).
    # For the additive group Z^2, this is calculated as g + h.
    result_deck = g + h

    print("Action 2: Deck Transformations (corresponds to Left Multiplication)")
    print("The action of g on h is calculated as g + h:")
    print(f"   ({g[0]}, {g[1]}) + ({h[0]}, {h[1]}) = ({result_deck[0]}, {result_deck[1]})")
    print("-" * 40)
    
    # Compare the results. The actions are the same if and only if the group is abelian.
    # Z^2 is an abelian group because vector addition is commutative.
    are_same = np.array_equal(result_holonomy, result_deck)
    
    print("Are the two actions the same for the torus T^2?")
    if are_same:
        print("Yes. The results are identical because the fundamental group of the torus, Z^2, is abelian (commutative).")
        print("For any g and h in Z^2, it is always true that h + g = g + h.")
    else:
        # This case will not be reached for the torus.
        print("No. This would imply the fundamental group is not abelian.")

# Run the demonstration.
demonstrate_actions_on_torus()