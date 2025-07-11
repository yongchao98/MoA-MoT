def calculate_riemann_components_kähler(m):
    """
    Calculates and prints the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension 'm'.

    The final formula is (m * (m + 1) / 2)^2.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension m must be a positive integer.")
        return

    # The real dimension n is 2 * m.
    n = 2 * m
    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {n}):")
    print("-" * 55)

    # --- Step 1: Calculate N ---
    # N represents the dimension of the space of symmetric (2,0)-tensors.
    # N = m * (m + 1) / 2
    n_val_numerator = m * (m + 1)
    n_val = int(n_val_numerator / 2)

    print(f"Step 1: Calculate N, the size of the equivalent Hermitian matrix.")
    print(f"Equation: N = (m * (m + 1)) / 2")
    print(f"Calculation: N = ({m} * ({m} + 1)) / 2 = {n_val_numerator} / 2 = {n_val}")

    # --- Step 2: Calculate N^2 ---
    # The number of independent real components in an N x N Hermitian matrix is N^2.
    result = n_val ** 2
    print(f"\nStep 2: Calculate the total number of independent components, which is N^2.")
    print(f"Equation: Components = N^2")
    print(f"Calculation: Components = {n_val}^2 = {result}")

    print("-" * 55)
    print(f"Result: The number of independent components is {result}.")
    print("\n")


if __name__ == '__main__':
    # Demonstrate the calculation for a few common dimensions.

    # For m=1 (e.g., a sphere, which is a Riemann surface).
    calculate_riemann_components_kähler(1)

    # For m=2 (e.g., a K3 surface or complex projective plane CP^2).
    calculate_riemann_components_kähler(2)

    # For m=3 (e.g., a typical Calabi-Yau 3-fold used in string theory).
    calculate_riemann_components_kähler(3)
