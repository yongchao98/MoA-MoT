import math

def display_nsvz_equation():
    """
    This function prints the NSVZ beta function equation and its components.
    The condition for this relation to hold is that the regularization scheme
    must preserve the holomorphy properties inherent to supersymmetry.
    """

    # --- Define symbolic representations for the equation's components ---
    beta_alpha = "β(α)"
    alpha = "α"
    T_G = "T_G"
    sum_Tr_gamma = "Σ_r T_r(1 - γ_r)"
    pi_symbol = "π"

    # --- Construct the equation string ---
    # The NSVZ beta function is an exact relation between the beta function
    # of the gauge coupling α = g^2 / (4π) and the anomalous dimensions γ_r
    # of the matter superfields.

    # Numerator contains the one-loop coefficient and the anomalous dimension contributions.
    # The number 3 comes from the gauge/gaugino loops.
    # The number 1 is from the (1 - gamma_r) term for each matter field.
    numerator = f"3*{T_G} - {sum_Tr_gamma}"

    # The denominator represents a resummation of certain higher-loop contributions.
    # The number 2 is part of the standard conversion from g to α.
    denominator = f"(1 - {T_G}*{alpha} / (2*{pi_symbol}))"

    # The coefficient in front sets the overall scale.
    # It contains the number 2.
    coefficient = f"-{alpha}**2 / (2*{pi_symbol})"

    # --- Print the full explanation and equation ---
    print("The exact condition for the NSVZ beta function to match non-renormalization theorems is that the chosen regularization scheme must preserve the holomorphy properties of the theory.")
    print("\nThis condition leads to the following exact equation:")
    print("-" * 50)
    print(f"{beta_alpha} = {coefficient} * [ {numerator} ] / [ {denominator} ]")
    print("-" * 50)

    print("\nWhere the terms represent:")
    print(f"  β(α): The beta function for the gauge coupling α.")
    print(f"  α: The fine-structure constant of the theory, α = g^2 / (4π).")
    print(f"  T_G: The Dynkin index for the adjoint representation (e.g., N for SU(N)). The coefficient '3' is explicit.")
    print(f"  Σ_r T_r: Sum of the Dynkin indices for all matter representations 'r'.")
    print(f"  γ_r: The anomalous dimension of the matter superfield in representation 'r'.")
    print(f"  π: The mathematical constant pi ({math.pi:.4f}...).")

display_nsvz_equation()