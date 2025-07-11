import sys

def demonstrate_torus_actions(gamma_mn, fiber_point_kl):
    """
    Demonstrates and compares the two actions of the fundamental group
    pi_1(T^2) on a fiber point of its universal cover.

    Args:
        gamma_mn (tuple): A tuple (m, n) representing an element of pi_1(T^2).
        fiber_point_kl (tuple): A tuple (k, l) representing a point in the fiber Z^2.
    """
    m, n = gamma_mn
    k, l = fiber_point_kl

    print(f"Let's consider the torus T^2.")
    print(f"The universal cover is R^2, and the fiber is the integer lattice Z^2.")
    print(f"We choose a loop gamma represented by the element ({m}, {n}) in pi_1(T^2).")
    print(f"We choose a point in the fiber p_inv(x_0) at ({k}, {l}).")
    print("-" * 40)

    # --- Action 1: Holonomy by lifting paths ---
    print("Action 1: Holonomy around loops (Path Lifting)")
    print(f"This action is defined by lifting the loop gamma starting from the fiber point ({k}, {l}).")
    
    # The lift of a path corresponding to (m, n) from (k, l) ends at (k+m, l+n).
    result1_k = k + m
    result1_l = l + n
    
    print(f"The calculation is a simple vector addition:")
    print(f"Result = ({k}, {l}) + ({m}, {n}) = ({k} + {m}, {l} + {n}) = ({result1_k}, {result1_l})")
    print("-" * 40)

    # --- Action 2: Restricting deck transformations ---
    print("Action 2: Restricting Deck Transformations to the Fiber")
    
    # The isomorphism between pi_1 and the Deck group is fixed by choosing a base point
    # in the fiber. We'll use (0, 0).
    fiber_base_point = (0, 0)
    print(f"1. First, we find the deck transformation g_gamma corresponding to the loop gamma=({m}, {n}).")
    print(f"   This is the transformation that maps the fiber base point {fiber_base_point} to the endpoint of the lift of gamma starting from {fiber_base_point}.")
    
    # The endpoint of the lift of gamma from (0, 0) is (m, n).
    # The deck transformations are translations g_(a,b)(x,y) = (x+a, y+b).
    # The transformation sending (0,0) to (m,n) is the translation by (m,n).
    deck_translation_vector = (m, n)
    print(f"   The lift of gamma from {fiber_base_point} ends at {deck_translation_vector}.")
    print(f"   Therefore, g_gamma is the translation by the vector ({m}, {n}).")
    
    print(f"\n2. Second, we apply this deck transformation to our fiber point ({k}, {l}).")
    
    # Apply the translation g_(m,n) to the point (k, l)
    result2_k = k + deck_translation_vector[0]
    result2_l = l + deck_translation_vector[1]

    print(f"   g_gamma(({k}, {l})) = ({k}, {l}) + ({m}, {n}) = ({k} + {m}, {l} + {n}) = ({result2_k}, {result2_l})")
    print("-" * 40)

    # --- Comparison ---
    print("Comparison:")
    print(f"Result of Action 1: ({result1_k}, {result1_l})")
    print(f"Result of Action 2: ({result2_k}, {result2_l})")
    
    if (result1_k, result1_l) == (result2_k, result2_l):
        print("\nThe two actions yield the same result.")
        print("This is because the fundamental group of the torus, Z^2, is abelian.")
    else:
        print("\nThe two actions yield different results.")

def main():
    # Example values
    # Let gamma be the loop corresponding to winding 2 times in the first direction
    # and 3 times in the second.
    gamma_element = (2, 3)

    # Let the point in the fiber be (5, 7).
    fiber_point = (5, 7)
    
    # You can change the input values here if you have Python installed.
    # For example:
    # if len(sys.argv) == 5:
    #   gamma_element = (int(sys.argv[1]), int(sys.argv[2]))
    #   fiber_point = (int(sys.argv[3]), int(sys.argv[4]))
    
    demonstrate_torus_actions(gamma_element, fiber_point)

if __name__ == "__main__":
    main()