import sys

def count_riemann_components_kahler(n):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given real dimension n.

    Args:
        n (int): The real dimension of the Kähler manifold.
    """
    print(f"Calculating the number of independent Riemann tensor components for a Kähler manifold of real dimension n = {n}.")

    # On a Kähler manifold, the real dimension n must be an even, positive integer.
    if not isinstance(n, int) or n <= 0 or n % 2 != 0:
        print("Error: The real dimension 'n' of a Kähler manifold must be a positive even integer.", file=sys.stderr)
        return

    # The complex dimension m is n / 2.
    m = n // 2
    print(f"The complex dimension is m = n / 2 = {n} / 2 = {m}.")
    print("-" * 30)

    # The formula for the number of components is (m * (m + 1) / 2)^2.
    print("The number of independent components is given by the formula:")
    print("N = (m * (m + 1) / 2)^2")
    print("\nPlugging in the value m =", m)

    # Step-by-step calculation
    step1_val = m + 1
    print(f"N = ({m} * ({m} + 1) / 2)^2")
    print(f"N = ({m} * {step1_val} / 2)^2")

    step2_val = m * step1_val
    print(f"N = ({step2_val} / 2)^2")

    d = step2_val // 2
    print(f"N = ({d})^2")

    result = d**2
    print(f"N = {result}")
    print("-" * 30)
    print(f"A Kähler manifold of real dimension {n} has {result} independent Riemann tensor components.")
    return result

# Example calculation for a 4-dimensional Kähler manifold (e.g., a K3 surface).
# The user can change this value to calculate for other dimensions.
real_dimension = 4
final_count = count_riemann_components_kahler(real_dimension)
# The final answer is wrapped according to the instruction format.
# This specific answer corresponds to the example case n=4.
# print(f"<<<{final_count}>>>")