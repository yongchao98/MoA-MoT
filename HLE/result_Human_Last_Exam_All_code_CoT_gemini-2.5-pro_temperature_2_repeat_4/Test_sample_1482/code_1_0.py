def formulate_it3_vertical_slice():
    """
    This function formulates and prints the mathematical equations for the
    vertical cross-section of a Gaussian Interval Type-3 Membership Function.
    """

    # --- Step 1: Define Symbolic Variables for Clarity ---
    # These represent the parameters of the Gaussian-based IT3 MF.
    x = "x"          # Primary input variable
    u = "u"          # Secondary input variable (for vertical slice, u in [0, 1])
    c = "c"          # Center of the Gaussian functions
    h = "h_L"        # Height of the lower slice of the Footprint of Uncertainty (FOU)
    s_U_upper = "σ_U_upper"  # Std. dev. for the Upper bound of the Upper FOU
    s_U_lower = "σ_U_lower"  # Std. dev. for the Lower bound of the Upper FOU
    s_L_upper = "σ_L_upper"  # Std. dev. for the Upper bound of the Lower FOU
    s_L_lower = "σ_L_lower"  # Std. dev. for the Lower bound of the Lower FOU

    # --- Step 2: Formulate the Four Bounding Gaussian MFs ---
    # These four functions define the boundaries of the entire IT3 fuzzy set.
    # We use descriptive names corresponding to standard literature.
    mu_UpperFOU_UpperBound = f"exp(-0.5 * (({x} - {c}) / {s_U_upper})^2)"
    mu_UpperFOU_LowerBound = f"exp(-0.5 * (({x} - {c}) / {s_U_lower})^2)"
    mu_LowerFOU_UpperBound = f"{h} * exp(-0.5 * (({x} - {c}) / {s_L_upper})^2)"
    mu_LowerFOU_LowerBound = f"{h} * exp(-0.5 * (({x} - {c}) / {s_L_lower})^2)"
    
    # --- Step 3: Formulate the Vertical Cross-Section Bounds ---
    # For a fixed 'x', the vertical cross-section is an IT2 MF over 'u'.
    # Its Upper and Lower Membership Functions are weighted averages of the
    # respective bounding MFs from Step 2, with 'u' as the weight.

    # Upper Membership Function of the vertical slice
    umf_slice_eq = f"{u} * ({mu_UpperFOU_UpperBound}) + (1 - {u}) * ({mu_UpperFOU_LowerBound})"

    # Lower Membership Function of the vertical slice
    lmf_slice_eq = f"{u} * ({mu_LowerFOU_UpperBound}) + (1 - {u}) * ({mu_LowerFOU_LowerBound})"

    # --- Step 4: Print the Final Formulation ---
    print("=" * 80)
    print("Formulation for the Vertical Cross-Section of a Gaussian IT3 Membership Function")
    print("=" * 80)
    print("\nAn Interval Type-3 Membership Function μ(x, u) is defined by its vertical")
    print("cross-sections. For a fixed primary input 'x', the cross-section is an")
    print("Interval Type-2 Membership Function defined over a secondary variable u ∈ [0, 1].")
    print("\nThis cross-section is characterized by its Upper Membership Function (UMF) and")
    print("Lower Membership Function (LMF), denoted as UMF_Slice(x, u) and LMF_Slice(x, u).")

    print("\n--- Foundational Bounding Functions (based on x) ---")
    print(f"\n1. Upper Bound of Upper FOU: μ_U_upper({x}) = {mu_UpperFOU_UpperBound}")
    print(f"\n2. Lower Bound of Upper FOU: μ_U_lower({x}) = {mu_UpperFOU_LowerBound}")
    print(f"\n3. Upper Bound of Lower FOU: μ_L_upper({x}) = {mu_LowerFOU_UpperBound}")
    print(f"\n4. Lower Bound of Lower FOU: μ_L_lower({x}) = {mu_LowerFOU_LowerBound}")
    print(f"\n   (Constraints: {s_U_upper} >= {s_U_lower}, {s_L_upper} >= {s_L_lower})")

    print("\n" + "=" * 80)
    print("--- Final Formulation for the Vertical Cross-Section (based on x and u) ---")

    print("\n1. Upper Membership Function of the Vertical Slice (UMF_Slice):\n")
    print(f"   UMF_Slice({x}, {u}) = {u} * μ_U_upper({x}) + (1 - {u}) * μ_U_lower({x})")
    print(f"   UMF_Slice({x}, {u}) = {umf_slice_eq}\n")

    print("\n2. Lower Membership Function of the Vertical Slice (LMF_Slice):\n")
    print(f"   LMF_Slice({x}, {u}) = {u} * μ_L_upper({x}) + (1 - {u}) * μ_L_lower({x})")
    print(f"   LMF_Slice({x}, {u}) = {lmf_slice_eq}\n")
    print("=" * 80)

# Execute the function to print the formulation.
if __name__ == "__main__":
    formulate_it3_vertical_slice()