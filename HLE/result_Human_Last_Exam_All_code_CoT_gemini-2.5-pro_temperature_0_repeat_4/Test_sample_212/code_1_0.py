def holonomy_action(gamma, fiber_point):
    """
    Calculates the action of a fundamental group element on a fiber point
    by holonomy (path lifting). For the torus T^2, this corresponds to
    vector addition.

    Args:
        gamma: A tuple (m, n) representing an element of pi_1(T^2).
        fiber_point: A tuple (m_prime, n_prime) representing a point in the fiber.

    Returns:
        A tuple representing the new fiber point.
    """
    m, n = gamma
    m_prime, n_prime = fiber_point
    
    # The action is defined by the endpoint of the lift of the loop gamma
    # starting at fiber_point. For the torus, this lift corresponds to adding
    # the vector representing the loop to the starting point vector.
    # The endpoint is (m_prime + m, n_prime + n).
    result = (m_prime + m, n_prime + n)
    return result

def deck_action(gamma, fiber_point):
    """
    Calculates the action of a fundamental group element on a fiber point
    by restricting the corresponding deck transformation.

    Args:
        gamma: A tuple (m, n) representing an element of pi_1(T^2).
        fiber_point: A tuple (m_prime, n_prime) representing a point in the fiber.

    Returns:
        A tuple representing the new fiber point.
    """
    m, n = gamma
    
    # 1. Find the deck transformation corresponding to gamma.
    # The base point in the fiber is chosen to be x_tilde_0 = (0, 0).
    # The deck transformation g_gamma maps x_tilde_0 to the endpoint of the
    # lift of gamma starting at x_tilde_0. This endpoint is (m, n).
    # For the torus, deck transformations are translations by integer vectors.
    # The transformation g_gamma is translation by the vector (m, n).
    translation_vector = gamma

    # 2. Apply this deck transformation (translation) to the fiber_point.
    m_prime, n_prime = fiber_point
    result = (m_prime + translation_vector[0], n_prime + translation_vector[1])
    return result

# --- Demonstration ---
# Let's choose an example element from the fundamental group and a point in the fiber.
gamma_element = (2, 3)
fiber_point_element = (5, 7)

# Calculate the result from the holonomy action
result_holonomy = holonomy_action(gamma_element, fiber_point_element)

# Calculate the result from the deck transformation action
result_deck = deck_action(gamma_element, fiber_point_element)

# --- Output ---
print(f"Let the fundamental group element be gamma = {gamma_element}")
print(f"Let the fiber point be x_tilde = {fiber_point_element}")
print("-" * 50)

print("Action 1: Holonomy around loops")
print(f"The new point is the result of the operation: {fiber_point_element} + {gamma_element}")
print(f"Result: ({fiber_point_element[0]} + {gamma_element[0]}, {fiber_point_element[1]} + {gamma_element[1]}) = {result_holonomy}")
print("-" * 50)

print("Action 2: Restricting deck transformations")
print(f"The deck transformation for gamma = {gamma_element} is translation by the vector {gamma_element}.")
print(f"Applying this to x_tilde = {fiber_point_element} gives: {fiber_point_element} + {gamma_element}")
print(f"Result: ({fiber_point_element[0]} + {gamma_element[0]}, {fiber_point_element[1]} + {gamma_element[1]}) = {result_deck}")
print("-" * 50)

# Final conclusion based on the demonstration
are_same = (result_holonomy == result_deck)
print(f"As the results are identical, the two actions are the same for T^2.")
print(f"Are the two actions the same? {are_same}")
