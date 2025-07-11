def calculate_riemann_components_kahler(m):
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of complex dimension m.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    # The formula for the number of components is (m * (m + 1) / 2)^2.
    # Let's calculate the intermediate part first.
    # This represents the dimension of the space of symmetric m x m matrices.
    d_val = m * (m + 1) // 2

    # The final number of components is the square of d_val.
    num_components = d_val ** 2

    print(f"For a Kähler manifold of complex dimension m = {m}:")
    print(f"The formula for the number of independent components is ( (m * (m + 1) / 2)^2 ).")
    print("\n--- Calculation Steps ---")
    print(f"1. First, we calculate the intermediate value 'd':")
    print(f"   d = ({m} * ({m} + 1)) / 2")
    print(f"   d = ({m * (m + 1)}) / 2")
    print(f"   d = {d_val}")
    print(f"\n2. Then, we square 'd' to get the final number of components:")
    print(f"   Number of components = d^2 = {d_val}^2")
    print(f"   Number of components = {num_components}")
    print("------------------------")


# --- Example Usage ---
# You can change this value to the complex dimension you are interested in.
complex_dimension = 3
calculate_riemann_components_kahler(complex_dimension)

# Another example for m=2
# calculate_riemann_components_kahler(2)