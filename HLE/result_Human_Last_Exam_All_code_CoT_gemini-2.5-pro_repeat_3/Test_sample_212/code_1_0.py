def holonomy_action(gamma, x_tilde):
    """
    Calculates the result of the holonomy action.
    
    Args:
        gamma (tuple): A pair of integers (m, n) representing an element of pi_1(T^2).
        x_tilde (tuple): A pair of integers (k, l) representing a point in the fiber.
        
    Returns:
        tuple: The resulting point in the fiber.
    """
    m, n = gamma
    k, l = x_tilde
    
    # The action is defined by the endpoint of the lifted path, which corresponds to vector addition.
    result_k = k + m
    result_l = l + n
    
    print("Action 1 (Holonomy):")
    print(f"Lifting the loop ({m}, {n}) starting from ({k}, {l}) ends at the point:")
    print(f"({k} + {m}, {l} + {n}) = ({result_k}, {result_l})")
    
    return (result_k, result_l)

def deck_action(gamma, x_tilde):
    """
    Calculates the result of the deck transformation action.

    Args:
        gamma (tuple): A pair of integers (m, n) representing an element of pi_1(T^2).
        x_tilde (tuple): A pair of integers (k, l) representing a point in the fiber.
        
    Returns:
        tuple: The resulting point in the fiber.
    """
    m, n = gamma
    k, l = x_tilde
    
    # The deck transformation corresponding to loop (m,n) is translation by vector (m,n).
    result_k = k + m
    result_l = l + n
    
    print("Action 2 (Deck Transformations):")
    print(f"Applying the deck transformation T_({m},{n}) to the point ({k}, {l}) results in:")
    print(f"({k} + {m}, {l} + {n}) = ({result_k}, {result_l})")
    
    return (result_k, result_l)

if __name__ == "__main__":
    # Choose an example element from the fundamental group pi_1(T^2)
    gamma_loop = (2, 3)
    
    # Choose an example point from the fiber p^-1(x_0)
    fiber_point = (5, 7)
    
    print(f"Let the loop gamma be represented by (m, n) = {gamma_loop}")
    print(f"Let the point in the fiber x_tilde be represented by (k, l) = {fiber_point}\n")
    
    result1 = holonomy_action(gamma_loop, fiber_point)
    print("-" * 20)
    result2 = deck_action(gamma_loop, fiber_point)
    print("-" * 20)
    
    if result1 == result2:
        print("\nThe results are identical. The two actions are the same for X = T^2.")
    else:
        print("\nThe results are different. The two actions are not the same for X = T^2.")
