import math

def solve_cylindrical_shell_problem():
    """
    Calculates the required permeability and internal field for a cylindrical
    magnetic shell to prevent external field distortion.

    This problem leads to a physically non-realizable solution for a non-trivial shell.
    The derivation shows that to prevent external field distortion, one of two
    conditions must be met:
    1. Trivial Case: The permeability of the shell is the same as the vacuum, mu = mu_0.
    2. Unphysical Case: The permeability is negative.

    The problem asks to exclude the trivial case, so we will calculate the result
    for the unphysical one, as this is the only other mathematical solution.
    """
    # Define the radii for the example calculation.
    # Let's use R1 = 1.0 and R2 = 2.0. The user can change these values.
    R1 = 1.0
    R2 = 2.0

    if R1 >= R2:
        print("Error: Inner radius R1 must be smaller than outer radius R2.")
        return

    R1_sq = R1**2
    R2_sq = R2**2

    # --- Part 1: Determine the required permeability (mu) ---
    # The condition for zero external distortion leads to the equation:
    # (mu^2 - mu_0^2) * (R2^2 - R1^2) = 0
    # For a non-trivial shell (R2 != R1), this gives mu^2 = mu_0^2.
    # The physical solution is mu = mu_0, which we must exclude.
    # The other mathematical solution arises from making the numerator of the
    # external field coefficient zero, which gives:
    # mu/mu_0 = - (R1^2 + R2^2) / (R2^2 - R1^2)
    # We calculate this unphysical relative permeability.

    mu_r_numerator = -(R1_sq + R2_sq)
    mu_r_denominator = (R2_sq - R1_sq)
    mu_r = mu_r_numerator / mu_r_denominator

    print("Required Permeability of the Shell Material (relative to vacuum permeability mu_0):")
    print("The only non-trivial solution requires a physically impossible negative permeability.")
    print("The formula for the relative permeability mu_r = mu/mu_0 is:")
    print(f"mu_r = - (R1^2 + R2^2) / (R2^2 - R1^2)")
    print("\nFor R1 = {} and R2 = {}:".format(R1, R2))
    print(f"mu/mu_0 = - ({R1_sq:.2f} + {R2_sq:.2f}) / ({R2_sq:.2f} - {R1_sq:.2f}) = {mu_r:.4f}")
    print("-" * 50)


    # --- Part 2: Find the magnetic field in the interior region (H_int) ---
    # With the value of mu found above, the internal field H_int is uniform
    # and aligned with the external field H0. Its magnitude is given by:
    # H_int / H0 = (R1^2 + R2^2) / R1^2

    h_int_ratio_numerator = (R1_sq + R2_sq)
    h_int_ratio_denominator = R1_sq
    h_int_ratio = h_int_ratio_numerator / h_int_ratio_denominator

    print("Magnetic Field in the Interior Region (H_int):")
    print("The field inside is uniform, given as a ratio to the applied field H0.")
    print("The formula for the ratio H_int / H0 is:")
    print(f"H_int / H0 = (R1^2 + R2^2) / R1^2")
    print("\nFor R1 = {} and R2 = {}:".format(R1, R2))
    print(f"H_int / H0 = ({R1_sq:.2f} + {R2_sq:.2f}) / {R1_sq:.2f} = {h_int_ratio:.4f}")
    print(f"So, H_int = {h_int_ratio:.4f} * H0 in the x-direction.")


# Execute the function to print the results.
solve_cylindrical_shell_problem()

# The permeability formula is mu/mu_0 = -(R1^2 + R2^2) / (R2^2 - R1^2)
# The internal field formula is H_int = H_0 * (R1^2 + R2^2) / R1^2
# For R1=1, R2=2: mu/mu_0 = -(1+4)/(4-1) = -5/3 = -1.6667
# H_int/H_0 = (1+4)/1 = 5
final_permeability = -(1**2 + 2**2) / (2**2 - 1**2)
final_H_ratio = (1**2 + 2**2) / (1**2)
answer = (final_permeability, final_H_ratio)