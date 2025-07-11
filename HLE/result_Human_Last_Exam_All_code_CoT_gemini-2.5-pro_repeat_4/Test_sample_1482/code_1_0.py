def generate_it3_mf_equation():
    """
    This function formulates and prints the equation for the upper bound
    of a vertical cross-section of a Gaussian Interval Type-3 Membership Function.
    """
    # --- Step 1: Define parameters for a specific example ---
    # These parameters define a specific Gaussian-based IT3 MF.
    # c: The center of the Gaussian function for the primary variable x.
    # sigma_L: The lower bound for the standard deviation range.
    # sigma_U: The upper bound for the standard deviation range.
    
    c = 5.0
    sigma_L = 1.0
    sigma_U = 2.0
    
    # --- Step 2: Explain the formulation step-by-step ---
    print("### Mathematical Formulation of an IT3 MF Vertical Cross-Section ###")
    print("\n1. Introduction:")
    print("An Interval Type-3 Membership Function (IT3 MF) models complex uncertainties. A vertical cross-section at a fixed primary input 'x' results in an Interval Type-2 (IT2) fuzzy set. This IT2 set is characterized by its own upper and lower membership functions, which depend on a secondary variable 'u' (where u is in [0, 1]).")

    print("\n2. Gaussian Paradigm:")
    print("We will use a Gaussian function as the base model. The general form is: exp(-0.5 * ((variable - center) / std_dev)^2).")
    
    print("\n3. Modeling Uncertainty in the Standard Deviation:")
    print("In an IT3 context, the parameters of the Gaussian, such as the standard deviation (σ), are functions of the secondary variable 'u'. For the upper bound of the membership function, we define an upper standard deviation function, σ_upper(u).")
    print("A common approach is to have σ_upper(u) linearly interpolate between a lower standard deviation (σ_L) and an upper one (σ_U) as 'u' varies from 0 to 1.")
    print(f"The formulation for this is: σ_upper(u) = u * σ_U + (1 - u) * σ_L")
    print(f"Using our parameters (σ_L={sigma_L}, σ_U={sigma_U}), this becomes:")
    print(f"σ_upper(u) = u * {sigma_U} + (1 - u) * {sigma_L}")
    
    print("\n4. Final Formulation:")
    print("By substituting the dynamic standard deviation σ_upper(u) into the Gaussian formula, we get the mathematical formulation for the upper bound of the vertical cross-section, denoted as μ_upper(x, u). This function defines the 'upper surface' of the domain of uncertainty.")

    # --- Step 3: Construct and print the final equation string ---
    # The final equation must have all numbers explicitly shown.
    final_equation_str = (
        f"μ_upper(x, u) = exp(-0.5 * ((x - {c}) / (u * {sigma_U} + (1.0 - u) * {sigma_L}))^2)"
    )
    
    print("\nThe final equation for the upper bound of the vertical cross-section is:")
    print(final_equation_str)

# Execute the function to display the explanation and the final equation.
generate_it3_mf_equation()