def count_riemann_tensor_components_kahler(m):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    # N is the dimension of the space of symmetric 2-tensors in m dimensions
    N_intermediate_product = m * (m + 1)
    N = N_intermediate_product / 2
    
    # The number of components is N^2
    num_components = N**2

    print(f"For a Kähler manifold of complex dimension m = {m}:")
    print(f"The number of independent components is given by the formula (m * (m + 1) / 2)^2.")
    # Show the calculation step-by-step
    print("Calculation:")
    print(f"({m} * ({m} + 1) / 2)^2")
    print(f"= ({m} * {m + 1} / 2)^2")
    print(f"= ({N_intermediate_product} / 2)^2")
    print(f"= ({int(N)})^2")
    print(f"= {int(num_components)}")

# As no dimension was specified, we use m=2 (a 4-dimensional real manifold) as an example.
complex_dimension_m = 2
count_riemann_tensor_components_kahler(complex_dimension_m)