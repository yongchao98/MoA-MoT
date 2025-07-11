def count_riemann_components_kahler(m: int):
    """
    Calculates and prints the number of independent real components of the
    Riemann tensor on a Kähler manifold of a given complex dimension.

    The formula for the number of components is (m * (m + 1) / 2)^2,
    where 'm' is the complex dimension of the manifold.

    Args:
        m: The complex dimension of the Kähler manifold (a positive integer).
    """
    if not isinstance(m, int) or m < 1:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    # The number of independent components is the square of the dimension of
    # the space of symmetric (2,0)-tensors.
    # The dimension of this space is N = m * (m + 1) / 2.
    # The number of independent real components of the Riemann tensor is N^2.

    m_plus_1 = m + 1
    numerator = m * m_plus_1
    N = numerator // 2
    num_components = N * N

    print(f"For a Kähler manifold of complex dimension m = {m}:")
    # The prompt requires printing each number in the final equation.
    print(f"The calculation is based on the formula: (m * (m + 1) / 2)^2")
    print(f"Step 1: m + 1 = {m} + 1 = {m_plus_1}")
    print(f"Step 2: m * (m + 1) = {m} * {m_plus_1} = {numerator}")
    print(f"Step 3: (m * (m + 1)) / 2 = {numerator} / 2 = {N}")
    print(f"Step 4: ({N})^2 = {N} * {N} = {num_components}")
    print(f"Result: There are {num_components} independent components.")
    print("-" * 30)

if __name__ == '__main__':
    print("This script calculates the number of independent components of the Riemann tensor")
    print("on a Kähler manifold of complex dimension 'm'.\n")

    # Demonstrate the calculation for a few example dimensions.
    # m=1 corresponds to a 2D real manifold like a sphere.
    count_riemann_components_kahler(1)

    # m=2 corresponds to a 4D real manifold like a K3 surface.
    count_riemann_components_kahler(2)

    # m=3 corresponds to a 6D real manifold like a Calabi-Yau 3-fold.
    count_riemann_components_kahler(3)