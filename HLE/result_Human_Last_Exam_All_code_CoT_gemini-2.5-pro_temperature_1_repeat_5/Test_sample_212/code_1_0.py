import numpy as np

def action_holonomy(gamma, fiber_point):
    """
    Calculates the action by holonomy (path lifting).
    For the torus, lifting a loop corresponding to vector gamma starting at fiber_point
    results in adding gamma to fiber_point.

    Args:
        gamma (tuple): A pair of integers (m, n) representing an element of pi_1(T^2).
        fiber_point (tuple): A pair of integers (k, l) representing a point in the fiber.

    Returns:
        tuple: The resulting point in the fiber.
    """
    m, n = gamma
    k, l = fiber_point
    result = (k + m, l + n)
    print("Action 1 (Holonomy):")
    print(f"Lifting the loop corresponding to {gamma} from the fiber point {fiber_point}...")
    print(f"The resulting endpoint is ({k} + {m}, {l} + {n}) = {result}\n")
    return result

def action_deck_transformation(gamma, fiber_point):
    """
    Calculates the action by restricting deck transformations.
    First, find the deck transformation for gamma, then apply it to fiber_point.

    Args:
        gamma (tuple): A pair of integers (m, n) representing an element of pi_1(T^2).
        fiber_point (tuple): A pair of integers (k, l) representing a point in the fiber.

    Returns:
        tuple: The resulting point in the fiber.
    """
    m, n = gamma
    k, l = fiber_point

    # The deck transformation for gamma = (m, n) is translation by (m, n).
    deck_transformation = lambda point: (point[0] + m, point[1] + n)
    
    result = deck_transformation(fiber_point)
    
    print("Action 2 (Deck Transformations):")
    print(f"The deck transformation for loop {gamma} is translation by vector {gamma}.")
    print(f"Applying this transformation to fiber point {fiber_point}...")
    print(f"The result is ({k} + {m}, {l} + {n}) = {result}\n")
    return result

if __name__ == "__main__":
    # Define an element of the fundamental group pi_1(T^2)
    # This corresponds to a loop wrapping 2 times horizontally and 3 times vertically.
    loop_gamma = (2, 3)

    # Define a point in the fiber p^-1(x_0)
    # This is a point on the integer lattice in R^2.
    fiber_point_x_tilde = (5, 7)

    print(f"Let X = T^2, the 2-torus.")
    print(f"The fundamental group pi_1(T^2) is Z^2.")
    print(f"The fiber p^-1(x_0) over the basepoint is the integer lattice Z^2.\n")
    print(f"Let's test the actions for a loop class {loop_gamma} and a fiber point {fiber_point_x_tilde}.")
    print("-" * 30)

    # Calculate the result of the first action
    result1 = action_holonomy(loop_gamma, fiber_point_x_tilde)

    # Calculate the result of the second action
    result2 = action_deck_transformation(loop_gamma, fiber_point_x_tilde)
    
    print("-" * 30)
    # Compare the results
    if result1 == result2:
        print("Conclusion: The results of the two actions are identical.")
    else:
        print("Conclusion: The results of the two actions are different.")
