import math

def calculate_min_zeros(chi_M, chi_dM):
    """
    Calculates the minimum number of zeros for a vector field on a compact
    manifold M with boundary dM.

    Args:
        chi_M (int): The Euler characteristic of the manifold M.
        chi_dM (int): The Euler characteristic of the boundary dM.
    """
    # The Euler characteristic of the boundary of a compact manifold is always even.
    if chi_dM % 2 != 0:
        print("Error: The Euler characteristic of a boundary (chi_dM) must be an even integer.")
        return

    # The formula for the minimum number of zeros is |χ(M) - χ(∂M)/2|
    # We use integer division // since chi_dM is guaranteed to be even.
    min_zeros = abs(chi_M - (chi_dM // 2))

    print("Given a compact manifold M with non-empty boundary ∂M:")
    print(f"  - Euler characteristic of M, χ(M) = {chi_M}")
    print(f"  - Euler characteristic of ∂M, χ(∂M) = {chi_dM}")
    print("\nThe minimum number of zeros a vector field on M can have is given by the formula: |χ(M) - χ(∂M)/2|")
    print("\nCalculation:")
    # The final code outputs each number in the final equation as requested.
    print(f"  |{chi_M} - {chi_dM}/2| = |{chi_M} - {chi_dM // 2}| = {min_zeros}")
    print(f"\nThe least number of zeros is {min_zeros}.")

if __name__ == '__main__':
    # Example: M is the 3-ball (D^3).
    # The Euler characteristic of D^3 is χ(M) = 1.
    # The boundary ∂M is the 2-sphere (S^2), with χ(∂M) = 2.
    print("--- Example: The 3-Ball (D^3) ---")
    chi_M_example = 1
    chi_dM_example = 2
    calculate_min_zeros(chi_M_example, chi_dM_example)

    # Example 2: M is the 2-disk (D^2).
    # The Euler characteristic of D^2 is χ(M) = 1.
    # The boundary ∂M is the 1-sphere (S^1), with χ(∂M) = 0.
    print("\n--- Example: The 2-Disk (D^2) ---")
    chi_M_example_2 = 1
    chi_dM_example_2 = 0
    calculate_min_zeros(chi_M_example_2, chi_dM_example_2)
