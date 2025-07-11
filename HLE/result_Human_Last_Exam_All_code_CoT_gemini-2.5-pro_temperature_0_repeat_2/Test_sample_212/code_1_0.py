def holonomy_action(gamma, point_in_fiber):
    """
    Calculates the action of a loop gamma on a point in the fiber
    using the holonomy (path lifting) definition.

    Args:
        gamma (tuple): A tuple (a, b) representing an element of pi_1(T^2).
        point_in_fiber (tuple): A tuple (m, n) representing a point in the fiber Z^2.

    Returns:
        tuple: The resulting point in the fiber.
    """
    a, b = gamma
    m, n = point_in_fiber
    
    # The lift of the path for gamma starting at (m, n) ends at (m+a, n+b).
    result_m = m + a
    result_n = n + b
    
    print("Holonomy Action:")
    print(f"  Loop gamma = {gamma}")
    print(f"  Starting point = {point_in_fiber}")
    print(f"  The lifted path ends at ({m} + {a}, {n} + {b}) = ({result_m}, {result_n})")
    
    return (result_m, result_n)

def deck_transformation_action(gamma, point_in_fiber, fiber_base_point=(0, 0)):
    """
    Calculates the action of a loop gamma on a point in the fiber
    by restricting the corresponding deck transformation.

    Args:
        gamma (tuple): A tuple (a, b) representing an element of pi_1(T^2).
        point_in_fiber (tuple): A tuple (m, n) representing a point in the fiber Z^2.
        fiber_base_point (tuple): The chosen base point in the fiber, usually (0,0).

    Returns:
        tuple: The resulting point in the fiber.
    """
    a, b = gamma
    m, n = point_in_fiber
    
    # Step 1: Find the deck transformation corresponding to gamma.
    # We lift gamma starting from the fiber_base_point.
    lift_endpoint_m = fiber_base_point[0] + a
    lift_endpoint_n = fiber_base_point[1] + b
    
    # The deck transformation is the translation that maps fiber_base_point to the endpoint.
    # deck_translation * fiber_base_point = lift_endpoint
    # So, the translation vector is (lift_endpoint - fiber_base_point), which is just gamma.
    deck_translation_vector = (a, b)
    
    # Step 2: Apply this deck transformation to the point_in_fiber.
    result_m = m + deck_translation_vector[0]
    result_n = n + deck_translation_vector[1]

    print("Deck Transformation Action:")
    print(f"  Loop gamma = {gamma}")
    print(f"  The corresponding deck transformation is a translation by {deck_translation_vector}.")
    print(f"  Applying this to the starting point {point_in_fiber}:")
    print(f"  Result is ({m} + {deck_translation_vector[0]}, {n} + {deck_translation_vector[1]}) = ({result_m}, {result_n})")

    return (result_m, result_n)

if __name__ == '__main__':
    # Define a sample loop in pi_1(T^2)
    # This corresponds to wrapping 3 times along the first circle and 5 times along the second.
    gamma_loop = (3, 5)
    
    # Define a sample point in the fiber p^-1(x_0) = Z^2
    fiber_point = (10, 20)
    
    print("--- Calculating Holonomy Action ---")
    result1 = holonomy_action(gamma_loop, fiber_point)
    print(f"Final point: {result1}\n")
    
    print("--- Calculating Deck Transformation Action ---")
    result2 = deck_transformation_action(gamma_loop, fiber_point)
    print(f"Final point: {result2}\n")
    
    # Compare the results
    print("--- Comparison ---")
    if result1 == result2:
        print("The results are the same.")
        print("Conclusion: For X = T^2, the two actions are the same.")
    else:
        print("The results are different.")
        print("Conclusion: For X = T^2, the two actions are not the same.")
