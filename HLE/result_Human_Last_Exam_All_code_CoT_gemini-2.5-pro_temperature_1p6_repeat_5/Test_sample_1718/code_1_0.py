import math

def calculate_riemann_components():
    """
    Calculates and prints the number of independent components of the Riemann tensor
    for a general Riemannian manifold and a Kähler manifold, given a complex dimension.
    """
    try:
        # Prompt user for the complex dimension
        m_str = input("Please enter the complex dimension (m) of the Kähler manifold: ")
        m = int(m_str)
        if m < 1:
            print("Error: The dimension must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter a whole number.")
        return

    # The real dimension is twice the complex dimension
    n = 2 * m

    print(f"\nFor a manifold of complex dimension m = {m}, the real dimension is n = {n}.")
    print("=====================================================================")

    # Case 1: General Riemannian Manifold
    print(f"1. For a general Riemannian manifold of real dimension n = {n}:")
    print(f"   The formula is: n^2 * (n^2 - 1) / 12")

    # Step-by-step calculation
    term1_gen = n**2
    term2_gen = term1_gen - 1
    numerator_gen = term1_gen * term2_gen
    result_gen = numerator_gen // 12

    print(f"   Calculation: ({n}^2 * ({n}^2 - 1)) / 12")
    print(f"              = ({term1_gen} * ({term1_gen} - 1)) / 12")
    print(f"              = ({term1_gen} * {term2_gen}) / 12")
    print(f"              = {numerator_gen} / 12")
    print(f"              = {result_gen}")
    print(f"   Result: There are {result_gen} independent components.")
    print("---------------------------------------------------------------------")

    # Case 2: Kähler Manifold
    print(f"2. For a Kähler manifold of complex dimension m = {m}:")
    print(f"   The formula is: ((m * (m + 1)) / 2)^2")

    # Step-by-step calculation
    term1_kahler = m + 1
    numerator_kahler = m * term1_kahler
    base_kahler = numerator_kahler // 2
    result_kahler = base_kahler**2
    
    print(f"   Calculation: (({m} * ({m} + 1)) / 2)^2")
    print(f"              = (({m} * {term1_kahler}) / 2)^2")
    print(f"              = ({numerator_kahler} / 2)^2")
    print(f"              = ({base_kahler})^2")
    print(f"              = {result_kahler}")
    print(f"   Result: There are {result_kahler} independent components.")
    print("=====================================================================")

if __name__ == "__main__":
    calculate_riemann_components()