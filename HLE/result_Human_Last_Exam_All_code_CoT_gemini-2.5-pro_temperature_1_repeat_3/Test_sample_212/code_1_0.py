def solve():
    """
    Demonstrates that for the torus T^2, the two actions of the fundamental group
    on the fiber of the universal cover are the same.
    """
    # Let's choose a point in the fiber p_inv(x_0). The fiber is the integer lattice Z^2.
    # A point is represented as a tuple (k, l).
    x_tilde = (5, 7)

    # Let's choose an element of the fundamental group pi_1(T^2).
    # The group is Z^2, and an element is represented as a tuple (m, n).
    gamma = (2, 3)

    print(f"Let's consider the torus T^2.")
    print(f"We choose a point in the fiber p_inv(x_0): x_tilde = {x_tilde}")
    print(f"We choose an element of the fundamental group pi_1(T^2): gamma = {gamma}\n")

    # --- Action 1: Holonomy (Path Lifting) ---
    print("1. Action by Holonomy (Path Lifting):")
    print("The action of gamma on x_tilde is found by lifting the loop gamma to a path")
    print("starting at x_tilde and finding its endpoint.")
    print("For the torus, this corresponds to vector addition on the universal cover R^2.")

    k, l = x_tilde
    m, n = gamma
    result1_k = k + m
    result1_l = l + n
    result1 = (result1_k, result1_l)

    print(f"New point = x_tilde + gamma")
    print(f"({k}, {l}) + ({m}, {n}) = ({k}+{m}, {l}+{n}) = {result1}")
    print(f"Result of Action 1: {result1}\n")

    # --- Action 2: Deck Transformations ---
    print("2. Action by Deck Transformations:")
    print(f"The deck transformation f_gamma corresponding to gamma = {gamma} is a translation by the vector {gamma}.")
    print(f"f_gamma(x, y) = (x + {m}, y + {n})")
    print(f"We apply this transformation to the point x_tilde = {x_tilde}.")
    
    # The deck transformation is defined by gamma
    def deck_transformation(point):
        return (point[0] + gamma[0], point[1] + gamma[1])

    result2 = deck_transformation(x_tilde)
    
    print(f"f_{gamma}({k}, {l}) = ({k}+{m}, {l}+{n}) = {result2}")
    print(f"Result of Action 2: {result2}\n")

    # --- Conclusion ---
    print("Conclusion:")
    if result1 == result2:
        print("The two actions produce the same result.")
        print("Therefore, for X = T^2, the actions are the same.")
    else:
        print("The two actions produce different results.")

solve()