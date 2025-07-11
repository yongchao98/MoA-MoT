import math

def calculate_riemann_components_kahler(m: int):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m: The complex dimension of the Kähler manifold.
    """
    if m <= 0:
        print(f"Dimension must be a positive integer, but got {m}.")
        return

    real_dimension = 2 * m
    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {real_dimension}):")

    # Step 1: Calculate the base of the square, D = m(m+1)/2
    d_numerator = m * (m + 1)
    
    # We use integer division as the numerator is always even
    d = d_numerator // 2

    # Step 2: Square the result
    n_components = d ** 2

    # Print the equation breakdown
    print(f"  The number of independent components is given by the formula (m * (m + 1) / 2)^2.")
    print(f"  Substituting m = {m}:")
    print(f"  N = ({m} * ({m} + 1) / 2)^2")
    print(f"    = ({m} * {m + 1} / 2)^2")
    print(f"    = ({d_numerator} / 2)^2")
    print(f"    = ({d})^2")
    print(f"    = {n_components}")
    print("-" * 20)

if __name__ == '__main__':
    # Calculate for the first four complex dimensions
    calculate_riemann_components_kahler(1)  # A Riemann surface
    calculate_riemann_components_kahler(2)  # A K3 surface or complex projective plane
    calculate_riemann_components_kahler(3)  # A Calabi-Yau threefold
    calculate_riemann_components_kahler(4)