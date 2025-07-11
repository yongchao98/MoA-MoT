import math

def solve_vector_field_zeros(chi_M, chi_dM):
    """
    Calculates the minimum number of zeros for a vector field on a compact manifold M
    with boundary dM, given their Euler characteristics.

    Args:
        chi_M (int): The Euler characteristic of the manifold M.
        chi_dM (int): The Euler characteristic of the boundary dM.
    """
    print(f"We are given a manifold M with chi(M) = {chi_M} and a boundary partial M with chi(partial M) = {chi_dM}.")
    print("The minimum number of zeros of a vector field on M can be determined by a condition relating these two values.")
    print("\nStep 1: Check if the manifold admits a non-vanishing vector field.")
    print("This is true if and only if 2 * chi(M) == chi(partial M).")
    
    # Check the condition
    if 2 * chi_M == chi_dM:
        result = 0
        print(f"\nStep 2: Evaluate the condition: 2 * {chi_M} == {chi_dM}")
        print(f"The condition is True, since {2 * chi_M} == {chi_dM}.")
        print("\nConclusion: The manifold admits a non-vanishing vector field.")
        print(f"The least number of zeros is 0.")
    else:
        result = abs(chi_M)
        print(f"\nStep 2: Evaluate the condition: 2 * {chi_M} == {chi_dM}")
        print(f"The condition is False, since {2 * chi_M} != {chi_dM}.")
        print("\nConclusion: The manifold does not admit a non-vanishing vector field.")
        print("In this case, the minimum number of zeros is given by the absolute value of chi(M).")
        print(f"The final equation is: |chi(M)| = |{chi_M}|")
        print(f"The least number of zeros is {result}.")

# Example: M is the 2-dimensional disk (D^2).
# chi(M) = 1
# The boundary, partial M, is a circle (S^1), so chi(partial M) = 0.
print("--- Example: M = 2-Disk (D^2) ---")
solve_vector_field_zeros(chi_M=1, chi_dM=0)

print("\n" + "="*40 + "\n")

# Example: M is the closed interval [0, 1].
# chi(M) = 1
# The boundary, partial M, is two points {0, 1}, so chi(partial M) = 2.
print("--- Example: M = Interval ([0, 1]) ---")
solve_vector_field_zeros(chi_M=1, chi_dM=2)