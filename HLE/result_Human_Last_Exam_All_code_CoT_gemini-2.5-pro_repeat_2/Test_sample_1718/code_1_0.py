import math

def calculate_kahler_riemann_components(m):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    n = 2 * m
    
    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {n}):")
    print("The number of independent real components (N) is given by the formula:")
    print("N = (m^2 * (m + 1)^2) / 4")
    print("\n--- Calculation Steps ---")

    # Step-by-step calculation
    m_plus_1 = m + 1
    m_sq = m**2
    m_plus_1_sq = m_plus_1**2
    numerator = m_sq * m_plus_1_sq
    denominator = 4
    result = numerator // denominator

    # Outputting the equation with numbers
    print(f"First, substitute m = {m} into the formula:")
    print(f"N = ({m}^2 * ({m} + 1)^2) / {denominator}")
    
    print(f"\nThen, evaluate the terms in the numerator:")
    print(f"N = ({m_sq} * {m_plus_1}^2) / {denominator}")
    print(f"N = ({m_sq} * {m_plus_1_sq}) / {denominator}")
    print(f"N = {numerator} / {denominator}")

    print(f"\nFinally, perform the division:")
    print(f"N = {result}")
    
    # This part is for the final answer block
    # print(f"\n<<<{result}>>>")


# We will use a 6-dimensional Kähler manifold (m=3) as an example.
calculate_kahler_riemann_components(3)