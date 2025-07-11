def holonomy_action(loop_element, fiber_point):
    """
    Calculates the action by holonomy (path lifting).
    For the torus T^2, a loop represented by the integer vector (a, b)
    acting on a fiber point (m, n) results in the endpoint of the lifted
    path, which is (m+a, n+b).
    """
    a, b = loop_element
    m, n = fiber_point
    result = (m + a, n + b)
    return result

def deck_transformation_action(loop_element, fiber_point):
    """
    Calculates the action by restricting deck transformations.
    For T^2, the loop (a, b) corresponds to the deck transformation
    phi(x, y) = (x+a, y+b). Applying this to the fiber point (m,n)
    gives the result.
    """
    a, b = loop_element
    m, n = fiber_point
    result = (m + a, n + b)
    return result

# Let's choose a specific loop from pi_1(T^2) and a point in the fiber.
# pi_1(T^2) is isomorphic to Z^2. A loop is represented by (a, b).
loop_gamma = (2, 3)
# The fiber p^-1(x_0) is isomorphic to Z^2. A fiber point is represented by (m, n).
fiber_point_x_tilde = (5, 7)

# Calculate the result of the first action (Holonomy)
holonomy_result = holonomy_action(loop_gamma, fiber_point_x_tilde)

# Calculate the result of the second action (Deck Transformations)
deck_result = deck_transformation_action(loop_gamma, fiber_point_x_tilde)

# Print the setup and the results of both actions to demonstrate they are identical.
print("We are examining two actions of the fundamental group on the fiber of the universal cover for the torus X = T^2.")
print(f"Let the fundamental group element be gamma = {loop_gamma}.")
print(f"Let the fiber point be x_tilde = {fiber_point_x_tilde}.")
print("-" * 30)

print("Action 1: By Holonomy (Path Lifting)")
print("The new point is the endpoint of the lift of gamma starting at x_tilde.")
a, b = loop_gamma
m, n = fiber_point_x_tilde
print(f"The calculation is adding the loop vector to the fiber point vector:")
print(f"Equation: ({m}, {n}) + ({a}, {b}) = ({m + a}, {n + b})")
print(f"Result: {holonomy_result}")
print("-" * 30)

print("Action 2: By Restricting Deck Transformations")
print("The loop gamma corresponds to a deck transformation phi_gamma, which we apply to x_tilde.")
print(f"The deck transformation for gamma=({a},{b}) is phi(x, y) = (x + {a}, y + {b}).")
print("Applying this to x_tilde:")
print(f"Equation: phi_({a},{b})({m}, {n}) = ({m} + {a}, {n} + {b})")
print(f"Result: {deck_result}")
print("-" * 30)

if holonomy_result == deck_result:
    print("Conclusion: The results are identical. The two actions are the same for T^2.")
else:
    print("Conclusion: The results are different. The two actions are not the same for T^2.")
