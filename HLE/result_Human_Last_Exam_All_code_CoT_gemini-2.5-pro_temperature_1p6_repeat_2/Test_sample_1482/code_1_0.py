def print_vertical_slice_formulation():
    """
    This function prints the mathematical formulation that characterizes the
    vertical cross-sections of a Gaussian Interval Type-3 Membership Function (IT3 MF).
    """

    print("The mathematical characterization of a vertical cross-section of an IT3 MF at a fixed primary input 'x' is defined by the Upper and Lower Membership Functions (UMF and LMF) of the resulting Interval Type-2 MF's Footprint of Uncertainty (FOU).")
    print("-" * 80)
    print("Variable Definitions:")
    print("  - x: The primary input variable.")
    print("  - u: The secondary input variable (i.e., the primary membership degree).")
    print("  - μ_upper_A_upper(x), μ_lower_A_upper(x): The UMF and LMF of the bounding upper IT2 MF.")
    print("  - μ_upper_A_lower(x), μ_lower_A_lower(x): The UMF and LMF of the bounding lower IT2 MF.")
    print("  - σ_u_upper, σ_u_lower: The standard deviations defining the shape of the Gaussian tails for the vertical slice's UMF and LMF, respectively.")
    print("-" * 80)

    # Formulation for the Upper Membership Function (UMF) of the Vertical Slice
    print("\nFormulation for the UMF of the Vertical Cross-Section, μ_upper_slice(u|x):")
    print("\nμ_upper_slice(u|x) = \n")
    
    # Piecewise definition
    print("  { 1,")
    print("  {    for u in the interval [ μ_lower_A_upper(x), μ_upper_A_upper(x) ]\n")
    
    print("  { exp( -( (u - μ_upper_A_upper(x)) ** 2 ) / (2 * (σ_u_upper ** 2)) ),")
    print("  {    for u > μ_upper_A_upper(x)\n")

    print("  { exp( -( (u - μ_lower_A_upper(x)) ** 2 ) / (2 * (σ_u_upper ** 2)) ),")
    print("  {    for u < μ_lower_A_upper(x)\n")

    print("-" * 80)

    # Formulation for the Lower Membership Function (LMF) of the Vertical Slice
    print("\nFormulation for the LMF of the Vertical Cross-Section, μ_lower_slice(u|x):")
    print("\nμ_lower_slice(u|x) = \n")

    # Piecewise definition
    print("  { 1,")
    print("  {    for u in the interval [ μ_lower_A_lower(x), μ_upper_A_lower(x) ]\n")

    print("  { exp( -( (u - μ_upper_A_lower(x)) ** 2 ) / (2 * (σ_u_lower ** 2)) ),")
    print("  {    for u > μ_upper_A_lower(x)\n")
    
    print("  { exp( -( (u - μ_lower_A_lower(x)) ** 2 ) / (2 * (σ_u_lower ** 2)) ),")
    print("  {    for u < μ_lower_A_lower(x)\n")


# Execute the function to print the formulation
print_vertical_slice_formulation()