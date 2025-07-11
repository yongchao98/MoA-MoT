import numpy as np

def action_by_deck_transformation(group_element, fiber_point):
    """
    Calculates the action of the fundamental group on the fiber
    by restricting the corresponding deck transformation.

    For the torus T^2, the universal cover is R^2, the fundamental group
    is Z^2, and the fiber over (0,0) is the integer lattice Z^2.

    A group element (m, n) from pi_1(T^2) corresponds to a deck
    transformation phi(x, y) = (x+m, y+n).

    Args:
        group_element (tuple): A tuple (m, n) representing an element of pi_1(T^2).
        fiber_point (tuple): A tuple (k, l) representing a point in the fiber p^-1(x_0).

    Returns:
        tuple: The resulting fiber point after the action.
    """
    m, n = group_element
    k, l = fiber_point
    
    # The action is vector addition
    result = (k + m, l + n)
    print("Action by Deck Transformation:")
    print(f"The deck transformation for group element {group_element} acts on fiber point {fiber_point}.")
    print(f"Result: {fiber_point} + {group_element} = {result}")
    return result

def action_by_holonomy(group_element, fiber_point):
    """
    Calculates the action of the fundamental group on the fiber
    by path lifting (holonomy/monodromy).

    A group element (m, n) corresponds to a loop on the torus.
    A simple representative is gamma(t) = (m*t mod 1, n*t mod 1).

    To find the action on a fiber point (k, l), we find the unique lift
    of this loop to a path in R^2 starting at (k, l). The lift is
    gamma_tilde(t) = (k + m*t, l + n*t).

    The action is the endpoint of this lifted path at t=1.

    Args:
        group_element (tuple): A tuple (m, n) representing an element of pi_1(T^2).
        fiber_point (tuple): A tuple (k, l) representing a point in the fiber p^-1(x_0).

    Returns:
        tuple: The resulting fiber point after the action.
    """
    m, n = group_element
    k, l = fiber_point
    
    # The action is the endpoint of the lifted path at t=1.
    t = 1
    result = (k + m * t, l + n * t)
    print("\nAction by Holonomy:")
    print(f"Lifting the loop for group element {group_element} starting from fiber point {fiber_point}.")
    print(f"The endpoint of the lift is ({k} + {m}*1, {l} + {n}*1) = {result}")
    return result

if __name__ == '__main__':
    # Choose a sample element from the fundamental group pi_1(T^2)
    g = (3, -2)

    # Choose a sample point from the fiber p^-1(x_0)
    y = (5, 8)

    print(f"Testing the two actions for X = T^2.")
    print(f"Fundamental group element g = {g}")
    print(f"Fiber point y = {y}")
    print("-" * 40)

    # Calculate the result of the first action
    result1 = action_by_deck_transformation(g, y)
    
    # Calculate the result of the second action
    result2 = action_by_holonomy(g, y)
    
    print("-" * 40)
    # Compare the results
    if result1 == result2:
        print("Conclusion: The two actions produce the same result.")
    else:
        print("Conclusion: The two actions produce different results.")
