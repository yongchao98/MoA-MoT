def formulate_it3_vertical_cross_section():
    """
    This function formulates and prints the mathematical expressions for the
    vertical cross-sections of an Interval Type-3 Membership Function (IT3 MF)
    using a Gaussian-based paradigm.
    """

    print("----------------------------------------------------------------------")
    print("Formulation for Vertical Cross-Sections of an IT3 MF")
    print("----------------------------------------------------------------------")
    print("\nAn Interval Type-3 Membership Function (IT3 MF) is denoted as μ̃(x, u).")
    print("For a fixed primary input x = x', the vertical cross-section μ̃(x', u)")
    print("is an Interval Type-2 Membership Function. This function is bounded by")
    print("an Upper Membership Function (UMF) and a Lower Membership Function (LMF).")
    print("\nUsing a Gaussian-based model, these bounds are formulated as follows:")

    # --- Variable Definitions ---
    print("\n--- Key Variables ---")
    print("μ̃(x', u): The vertical cross-section at a fixed primary input x'.")
    print("μ̄_μ̃(x', u): The Upper Membership Function (UMF) of the vertical cross-section.")
    print("μ̲_μ̃(x', u): The Lower Membership Function (LMF) of the vertical cross-section.")
    print("u          : The secondary input variable (primary membership grade), where u ∈ [0, 1].")
    print("c_μ̄(x')    : The center (mean) of the Gaussian for the UMF at x'.")
    print("σ_μ̄(x')    : The standard deviation of the Gaussian for the UMF at x'.")
    print("c_μ̲(x')    : The center (mean) of the Gaussian for the LMF at x'.")
    print("σ_μ̲(x')    : The standard deviation of the Gaussian for the LMF at x'.")

    # --- Equation for the Upper Membership Function (UMF) ---
    print("\n--- Equation for the Upper Membership Function (UMF) ---")
    # Using unicode characters for better readability
    equation_umf = "μ̄_μ̃(x', u) = exp( - (u - c_μ̄(x'))² / (2 * σ_μ̄(x')²) )"
    print(equation_umf)
    print("\nBreaking down the UMF equation:")
    print(f"  - Expression: exp(...) represents Euler's number 'e' raised to the power of the expression.")
    print(f"  - Numerator:  -(u - c_μ̄(x'))² is the squared difference between the secondary variable 'u' and the UMF center.")
    print(f"  - Denominator: (2 * σ_μ̄(x')²) is twice the squared standard deviation of the UMF, controlling the spread.")

    # --- Equation for the Lower Membership Function (LMF) ---
    print("\n--- Equation for the Lower Membership Function (LMF) ---")
    equation_lmf = "μ̲_μ̃(x', u) = exp( - (u - c_μ̲(x'))² / (2 * σ_μ̲(x')²) )"
    print(equation_lmf)
    print("\nBreaking down the LMF equation:")
    print(f"  - Expression: exp(...) represents Euler's number 'e' raised to the power of the expression.")
    print(f"  - Numerator:  -(u - c_μ̲(x'))² is the squared difference between the secondary variable 'u' and the LMF center.")
    print(f"  - Denominator: (2 * σ_μ̲(x')²) is twice the squared standard deviation of the LMF, controlling its spread.")
    print("----------------------------------------------------------------------")


if __name__ == '__main__':
    formulate_it3_vertical_cross_section()