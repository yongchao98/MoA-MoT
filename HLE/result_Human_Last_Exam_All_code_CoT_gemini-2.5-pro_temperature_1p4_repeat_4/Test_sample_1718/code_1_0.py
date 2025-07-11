def calculate_riemann_components(m):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Please provide a positive integer for the complex dimension.")
        return

    print(f"Calculating for a Kähler manifold of complex dimension m = {m}")
    print("----------------------------------------------------------")

    # The formula for the number of components is ( (m * (m + 1)) / 2 )^2.
    # We break it down to show the intermediate steps.

    # Step 1: Calculate N, the dimension of the space of symmetric 2-tensors.
    # This corresponds to the size of the equivalent Hermitian matrix.
    m_plus_1 = m + 1
    numerator_N = m * m_plus_1
    # We use integer division // as the result is always an integer.
    N = numerator_N // 2

    print(f"The calculation relies on the formula: ( (m * (m + 1)) / 2 )^2")
    print(f"Step 1: Calculate the term in the parentheses, let's call it N.")
    print(f"   N = ({m} * ({m} + 1)) / 2")
    print(f"   N = ({m} * {m_plus_1}) / 2")
    print(f"   N = {numerator_N} / 2 = {N}")

    # Step 2: The number of independent components is N^2.
    result = N * N
    print(f"Step 2: The number of independent components is N^2.")
    print(f"   Number of components = {N}^2 = {result}")
    print(f"Result: The number of independent entries is {result}.\n")

if __name__ == '__main__':
    # Demonstrate the calculation for a few example dimensions.
    # m=1 corresponds to a Riemann surface (real 2-manifold).
    calculate_riemann_components(1)

    # m=2 corresponds to a K3 surface or other complex surface (real 4-manifold).
    calculate_riemann_components(2)

    # m=3 corresponds to a Calabi-Yau manifold (real 6-manifold).
    calculate_riemann_components(3)