import math

def calculate_nsvz_beta_sqed(alpha, Nf):
    """
    Calculates the NSVZ beta function for Supersymmetric QED (SQED).

    In SQED (a U(1) gauge theory), the gauge group is abelian, so the
    adjoint Casimir T(adj) is 0. We consider Nf pairs of chiral superfields
    with charge +/- 1.

    The NSVZ beta function formula is:
    β(α) = - (α² / 2π) * [3T(adj) - Σ T(R_i)(1 - γ_i)] / [1 - T(adj) * α / π]

    For SQED (U(1) gauge group, Nf flavors):
    - T(adj) = 0
    - Σ T(R_i) = 2 * Nf (sum over Nf fields and Nf anti-fields)
    - The anomalous dimension for each field is γ_i = -α/π.

    This simplifies the formula to:
    β(α) = - (α² / 2π) * [-2 * Nf * (1 - γ)] = (Nf * α² / π) * (1 - γ)
    """

    # The anomalous dimension in the NSVZ scheme for a field in SQED
    gamma = -alpha / math.pi

    # The NSVZ beta function for SQED
    beta = (Nf * alpha**2 / math.pi) * (1 - gamma)

    print(f"Calculating the NSVZ beta function for SQED with parameters:")
    print(f"Fine-structure constant alpha = {alpha}")
    print(f"Number of flavors Nf = {Nf}")
    print("-" * 20)

    # Output each number in the final equation as requested
    print("Step 1: Calculate the anomalous dimension (γ)")
    print(f"γ = -α / π = -{alpha} / {math.pi:.5f} = {gamma:.5f}")
    print("")

    print("Step 2: Substitute into the NSVZ formula for SQED")
    print("β(α) = (Nf * α² / π) * (1 - γ)")
    print(f"β({alpha}) = ({Nf} * {alpha}² / {math.pi:.5f}) * (1 - ({gamma:.5f}))")

    # Breaking down the calculation
    term1_val = Nf * alpha**2 / math.pi
    term2_val = 1 - gamma
    print(f"β({alpha}) = ({term1_val:.5f}) * ({term2_val:.5f})")
    print("-" * 20)
    print(f"Final Result: β({alpha}) = {beta:.5f}")


# --- User Parameters ---
# You can change these values to see how the beta function changes.
fine_structure_constant = 0.1
number_of_flavors = 3
# --- End of Parameters ---

calculate_nsvz_beta_sqed(fine_structure_constant, number_of_flavors)
