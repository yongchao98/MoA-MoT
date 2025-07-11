def get_least_number_of_zeros(chi_M):
    """
    Calculates the least number of zeros a vector field can have on a compact
    manifold M with Euler characteristic chi_M.

    The result is based on a theorem by Hopf: a nowhere-zero vector field
    exists on a compact manifold if and only if its Euler characteristic is 0.
    If the Euler characteristic is non-zero, a vector field with a single
    zero can always be constructed.

    Args:
        chi_M (int): The Euler characteristic of the manifold M.

    Returns:
        int: The minimum number of zeros.
    """
    if chi_M == 0:
        # If chi(M) is 0, a nowhere-zero vector field exists.
        min_zeros = 0
    else:
        # If chi(M) is not 0, any vector field must have at least one zero,
        # and a field with exactly one zero can be constructed.
        min_zeros = 1
    return min_zeros

# --- Examples ---

# Example 1: A cylinder or a Mobius strip.
# For a cylinder, chi(M) = 0.
chi_M_cylinder = 0
min_zeros_cylinder = get_least_number_of_zeros(chi_M_cylinder)
print(f"For a manifold M with chi(M) = {chi_M_cylinder}, the least number of zeros a vector field can have is:")
print(min_zeros_cylinder)

print("-" * 20)

# Example 2: A disk (or any contractible manifold with boundary).
# For a disk, chi(M) = 1.
chi_M_disk = 1
min_zeros_disk = get_least_number_of_zeros(chi_M_disk)
print(f"For a manifold M with chi(M) = {chi_M_disk}, the least number of zeros a vector field can have is:")
print(min_zeros_disk)

print("-" * 20)

# Example 3: A sphere with 4 holes (a "four-holed sphere").
# chi(M) = 2 - 2*g - b = 2 - 2*0 - 4 = -2.
chi_M_surface = -2
min_zeros_surface = get_least_number_of_zeros(chi_M_surface)
print(f"For a manifold M with chi(M) = {chi_M_surface}, the least number of zeros a vector field can have is:")
print(min_zeros_surface)
