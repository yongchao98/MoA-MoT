def count_riemann_tensor_components_on_kahler(complex_dim):
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        complex_dim (int): The complex dimension 'm' of the Kähler manifold.
                           Must be a non-negative integer.
    """
    if not isinstance(complex_dim, int) or complex_dim < 0:
        print("Error: The complex dimension must be a non-negative integer.")
        return

    m = complex_dim

    # The formula for the number of components is (m * (m + 1) / 2)^2.
    # First, we calculate the base of the square, N = m * (m+1) / 2
    numerator = m * (m + 1)
    
    # Since either m or m+1 is even, the numerator is always even.
    # We can safely use integer division.
    base = numerator // 2
    
    # The total number of components is N^2.
    num_components = base ** 2

    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {2*m}):")
    print(f"The number of independent components is given by the formula N = (m * (m + 1) / 2)^2.")
    print(f"The calculation proceeds as follows:")
    # We explicitly show each number in the equation
    print(f"N = ({m} * ({m} + 1) / 2)^2")
    print(f"N = ({m} * {m+1} / 2)^2")
    print(f"N = ({numerator} / 2)^2")
    print(f"N = {base}^2")
    print(f"N = {num_components}")
    print("-" * 30)

# --- Example Calculations ---

# For m=1 (a 2-dimensional real manifold, like a sphere).
# Any 2D Riemannian manifold is Kähler. The number of components should be 1.
count_riemann_tensor_components_on_kahler(1)

# For m=2 (a 4-dimensional real manifold, like a K3 surface).
# A general 4D Riemannian manifold has 20 components. A Kähler 4-manifold has fewer.
count_riemann_tensor_components_on_kahler(2)

# For m=3 (a 6-dimensional real manifold, like a Calabi-Yau 3-fold).
count_riemann_tensor_components_on_kahler(3)