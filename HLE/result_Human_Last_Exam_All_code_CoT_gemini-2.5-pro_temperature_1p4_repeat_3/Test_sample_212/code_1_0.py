import sys

def solve():
    """
    This script demonstrates that for X = T^2, the two actions of the
    fundamental group on the fiber of the universal cover are the same.
    """
    # Define parameters for the demonstration.
    # An element of the fundamental group pi_1(T^2) ~= Z^2. Let's use g = (2, 1).
    m, n = 2, 1
    # A point in the fiber p^-1(x_0) ~= Z^2. Let's use x_tilde = (3, 4).
    k, l = 3, 4

    # Represent the group element and fiber points as tuples.
    g = (m, n)
    x_tilde = (k, l)
    # The standard base point in the fiber is (0, 0).
    x_tilde_0 = (0, 0)

    print("We will verify that the two actions of the fundamental group on the fiber are the same for the torus T^2.")
    print("The fundamental group pi_1(T^2) is isomorphic to Z^2.")
    print("The fiber p^-1(x_0) is the integer lattice Z^2.")
    print("-" * 60)
    print(f"Let's choose a group element g = {g} from pi_1(T^2).")
    print(f"Let's choose a fiber point x_tilde = {x_tilde} from the fiber.")
    print(f"The chosen base point in the fiber is x_tilde_0 = {x_tilde_0}.")
    print("-" * 60)

    # --- Action 1: Holonomy Action ---
    print("Action 1: The action by holonomy (path lifting).")
    print("This action lifts the loop for g to a path starting at x_tilde. The result is the endpoint of this path.")
    print(f"For the torus, this corresponds to adding the vector for g to the vector for x_tilde.")
    
    result1_k = k + m
    result1_l = l + n
    result1 = (result1_k, result1_l)
    
    print(f"Final Equation (Action 1): {x_tilde} + {g} = {result1}")
    print(f"Result of Action 1: {result1}")
    print("-" * 60)

    # --- Action 2: Deck Transformation Action ---
    print("Action 2: The action by restricting deck transformations.")
    print("This involves two steps:")
    print("  a) Identify the group element g with a deck transformation phi_g.")
    print("  b) Apply phi_g to the fiber point x_tilde.")
    
    print("\nStep a: Find the deck transformation phi_g for g.")
    print(f"phi_g is the unique deck transformation that maps the fiber base point {x_tilde_0} to the endpoint of the path for g lifted from {x_tilde_0}.")
    
    lift_endpoint_k = x_tilde_0[0] + g[0]
    lift_endpoint_l = x_tilde_0[1] + g[1]
    lift_endpoint = (lift_endpoint_k, lift_endpoint_l)
    print(f"The endpoint of the lifted path is {x_tilde_0} + {g} = {lift_endpoint}.")
    
    print(f"Deck transformations for the torus are translations. So, phi_g must be the translation by the vector {lift_endpoint}.")
    
    print("\nStep b: Apply phi_g to x_tilde.")
    print(f"phi_g acts on x_tilde by adding the translation vector.")
    
    result2_k = x_tilde[0] + lift_endpoint[0]
    result2_l = x_tilde[1] + lift_endpoint[1]
    result2 = (result2_k, result2_l)

    print(f"Final Equation (Action 2): {x_tilde} + {lift_endpoint} = {result2}")
    print(f"Result of Action 2: {result2}")
    print("-" * 60)

    # --- Conclusion ---
    print("Comparison:")
    print(f"Action 1 Result: {result1}")
    print(f"Action 2 Result: {result2}")

    if result1 == result2:
        print("\nThe two actions produce the same result.")
        print("This is true in general for the torus because its fundamental group, Z^2, is abelian.")
    else:
        # This case should not be reached for the torus.
        print("\nThe two actions produce different results.")

solve()