import numpy as np

def holonomy_action(group_element, fiber_point):
    """
    Computes the action of a fundamental group element on a fiber point
    via holonomy (path lifting).
    
    For the torus T^2, this corresponds to vector addition.
    
    Args:
        group_element (tuple of int): A pair (a, b) representing an element of pi_1(T^2).
        fiber_point (tuple of int): A pair (m, n) representing a point in the fiber Z^2.
        
    Returns:
        tuple of int: The new fiber point after the action.
    """
    # The endpoint of the lifted path starting at fiber_point is fiber_point + group_element
    result = np.add(fiber_point, group_element)
    return tuple(result)

def deck_action(group_element, fiber_point):
    """
    Computes the action of a fundamental group element on a fiber point
    via the corresponding deck transformation.
    
    Args:
        group_element (tuple of int): A pair (a, b) representing an element of pi_1(T^2).
        fiber_point (tuple of int): A pair (m, n) representing a point in the fiber Z^2.
        
    Returns:
        tuple of int: The new fiber point after the action.
    """
    # The deck transformation corresponding to group_element (a, b) is a translation by that vector.
    # Applying this translation to the fiber_point gives the result.
    deck_transformation = lambda point: np.add(point, group_element)
    result = deck_transformation(fiber_point)
    return tuple(result)

def main():
    # Let's choose a sample point in the fiber p^{-1}(x_0)
    fiber_point = (5, 2)
    
    # And a sample element from the fundamental group pi_1(T^2)
    group_element = (3, 4)
    
    print("Let's test the two actions for X = T^2.\n")
    print(f"Fiber point f = {fiber_point}")
    print(f"Fundamental group element g = {group_element}\n")

    # --- Calculate Action 1: Holonomy ---
    result_holonomy = holonomy_action(group_element, fiber_point)
    print("Action 1: Holonomy")
    # Output each number in the final equation
    print(f"f' = f + g = ({fiber_point[0]}, {fiber_point[1]}) + ({group_element[0]}, {group_element[1]}) = ({result_holonomy[0]}, {result_holonomy[1]})")

    print("-" * 20)

    # --- Calculate Action 2: Deck Transformation ---
    result_deck = deck_action(group_element, fiber_point)
    print("Action 2: Deck Transformation")
    # Output each number in the final equation
    print(f"The deck transformation for g is phi_g(p) = p + g.")
    print(f"f' = phi_g(f) = ({fiber_point[0]}, {fiber_point[1]}) + ({group_element[0]}, {group_element[1]}) = ({result_deck[0]}, {result_deck[1]})")
    
    print("\n" + "=" * 20)
    # --- Compare Results ---
    are_same = (result_holonomy == result_deck)
    print(f"Are the resulting fiber points the same? {are_same}")
    print("Conclusion: For the 2-torus, the two actions are indeed the same.")

if __name__ == "__main__":
    main()