import sys

def solve():
    """
    This script demonstrates that for the 2-torus T^2, the two descriptions
    of the action of the fundamental group pi_1(T^2) on the fiber of its
    universal cover are indeed the same.

    We use the following setup:
    - The torus X = T^2 is represented as R^2 / Z^2.
    - The universal cover is X_tilde = R^2.
    - The covering map is p(x,y) = (x mod 1, y mod 1).
    - The base point in X is x0 = p(0,0).
    - The fiber over x0 is the integer lattice p^-1(x0) = Z^2.
    - The fundamental group pi_1(T^2, x0) is isomorphic to Z^2.
    """

    # Let's choose a point in the fiber p^-1(x0). A point in the fiber is a pair of integers.
    # We can represent it as a tuple (a, b).
    fiber_point = (3, 5)

    # Let's choose an element of the fundamental group pi_1(T^2, x0).
    # This is also represented by a pair of integers (m, n).
    # A representative loop is gamma(t) = (m*t mod 1, n*t mod 1) for t in [0,1].
    loop_class = (2, 1)

    print(f"We are considering the case where X is the torus T^2.")
    print(f"The universal cover is R^2, and the fiber over the basepoint is the integer lattice Z^2.")
    print(f"We will test the two actions on the fiber point x_tilde = {fiber_point}.")
    print(f"The element of the fundamental group is represented by the vector {loop_class}.\n")

    a, b = fiber_point
    m, n = loop_class

    # --- Action 1: Holonomy around a loop ---
    print("--- Action 1: Holonomy (Monodromy) Action ---")
    print("This action is defined by lifting the loop to the universal cover and finding the endpoint.")

    # The lift of the loop gamma corresponding to (m,n) starting at (a,b) is
    # the path gamma_tilde(t) = (a + m*t, b + n*t) in R^2.
    # The action of the loop class on the fiber point is the endpoint of this path at t=1.
    endpoint_a = a + m
    endpoint_b = b + n
    result_action1 = (endpoint_a, endpoint_b)

    print(f"The loop corresponding to {loop_class} is lifted to a path in R^2 starting at {fiber_point}.")
    print(f"The lifted path is gamma_tilde(t) = ({a} + {m}*t, {b} + {n}*t).")
    print(f"The endpoint of this path at t=1 defines the action.")
    print(f"Resulting point: ({a} + {m}, {b} + {n}) = {result_action1}.")
    print(f"Result of Action 1: {result_action1}\n")


    # --- Action 2: Restricting deck transformations ---
    print("--- Action 2: Deck Transformation Action ---")
    print("This action is defined by finding the deck transformation corresponding to the loop, and then applying it to the fiber point.")

    # Step 2a: Find the deck transformation corresponding to the loop class (m,n).
    # We do this by fixing a basepoint in the fiber, typically (0,0), let's call it x_tilde_0.
    fiber_basepoint = (0, 0)

    # The lift of the loop gamma from x_tilde_0 is (m*t, n*t). The endpoint is (m, n) at t=1.
    endpoint_from_origin = (m, n)

    # For the cover R^2 -> T^2, deck transformations are translations by integer vectors.
    # The transformation phi(x,y) that sends (0,0) to (m,n) is phi(x,y) = (x+m, y+n).
    # So the deck transformation is a translation by the vector (m,n), which represents the loop class.
    deck_translation_vector = loop_class
    print(f"The deck transformation for the loop class {loop_class} is found by seeing where it sends a basepoint of the fiber, e.g., (0,0).")
    print(f"The lift from (0,0) ends at {endpoint_from_origin}, so the deck transformation is a translation by the vector {deck_translation_vector}.")

    # Step 2b: Apply this deck transformation to our chosen fiber point.
    result_action2_a = a + m
    result_action2_b = b + n
    result_action2 = (result_action2_a, result_action2_b)

    print(f"Applying this translation to our fiber point {fiber_point}:")
    print(f"Resulting point: ({a} + {m}, {b} + {n}) = {result_action2}.")
    print(f"Result of Action 2: {result_action2}\n")


    # --- Conclusion ---
    print("--- Conclusion ---")
    if result_action1 == result_action2:
        print("The results of the two actions are identical.")
        answer = "Yes"
    else:
        # This case will not be reached for universal covers.
        print("The results of the two actions are different.")
        answer = "No"

    print(f"\nAre these two actions the same when X = T^2? {answer}.")

if __name__ == "__main__":
    solve()