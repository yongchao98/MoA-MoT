def explain_ideal_ni_ce_ratio():
    """
    Explains the factors influencing the ideal Ni/Ce ratio for catalysis
    and provides a literature-based range.
    """

    print("Analyzing the ideal Ni/Ce ratio for Ni-Ceria catalysts...")
    print("-" * 60)

    # --- Step 1: Explain the determining factors ---
    print("Factor 1: The ideal ratio is not a single value. It depends on:")
    print("  - Synthesis method (e.g., co-precipitation, impregnation)")
    print("  - Operating temperature and pressure")
    print("  - Specific reaction (Water Gas Shift vs. Water Splitting)")
    print("\n")

    # --- Step 2: Explain the underlying scientific principle ---
    print("Principle: Catalytic performance is maximized by increasing the number of active sites.")
    print("This is achieved with high dispersion of Ni on the Ceria support, creating a strong Ni-O-Ce interface.")
    print("High Ni content can lead to particle agglomeration (sintering), which reduces catalytic activity.")
    print("\n")

    # --- Step 3: Present the literature-supported optimal range ---
    print("Conclusion from Literature Review:")
    print("For both Water Gas Shift (WGS) and Water Splitting (WS) reactions, a relatively low Ni content is generally preferred.")
    
    # Define the range for the final equation
    lower_bound_atomic_percent = 5
    upper_bound_atomic_percent = 15
    
    lower_bound_ratio = 0.05
    upper_bound_ratio = 0.15
    
    print(f"Optimal performance is often found when the Ni atomic percentage is between {lower_bound_atomic_percent}% and {upper_bound_atomic_percent}%.")
    
    # Final equation format
    print("\nThis corresponds to the following equation for the atomic ratio R = Ni / (Ni + Ce):")
    print(f"{lower_bound_ratio} <= R <= {upper_bound_ratio}")
    print("-" * 60)


# Run the explanatory function
explain_ideal_ni_ce_ratio()
