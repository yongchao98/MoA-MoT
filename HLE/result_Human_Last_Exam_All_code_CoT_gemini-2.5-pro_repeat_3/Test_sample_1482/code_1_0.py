import math

def generate_it3_formulation():
    """
    Generates and prints the mathematical formulation for the vertical cross-section
    of an Interval Type-3 Membership Function using a Gaussian paradigm.
    """

    # --- Step 1: Define parameters for the Gaussian functions ---
    # These are example values to make the final equation concrete.
    # For the primary Upper Membership Function (UMF)
    c_upper = 5.0  # Mean of the primary Gaussian
    sigma_upper = 1.0  # Standard deviation of the primary Gaussian

    # For the vertical cross-section's shape
    beta = 0.2  # Proportionality constant for the spread

    # --- Step 2: Construct the mathematical expressions as strings ---

    # Expression for the primary UMF, which determines the center of the vertical slice
    primary_umf_str = f"exp(-0.5 * ((x - {c_upper}) / {sigma_upper})^2)"

    # Expression for the upper bound of the vertical cross-section at a fixed x.
    # This is a function of the secondary variable 'u'.
    # The notation μ_upper(A(x))(u) represents the upper membership grade of u
    # in the IT2 MF that is the vertical slice of the IT3 MF A at x.
    vertical_slice_umf_str = (
        f"μ_upper(A(x))(u) = exp(-0.5 * ((u - μ_upper_A(x)) / "
        f"({beta} * μ_upper_A(x)))^2)"
    )

    # --- Step 3: Print the final formulation ---
    print("The mathematical formulation for the upper bound of a vertical cross-section of an IT3 MF is given in two parts:")
    print("\nPart 1: The Primary Upper Membership Function (UMF), which defines the domain of the secondary variable.")
    print(f"μ_upper_A(x) = {primary_umf_str}")

    print("\nPart 2: The Upper Membership Function of the Vertical Cross-Section at a fixed 'x'.")
    print("This function defines the shape of the uncertainty in the secondary dimension 'u'.")
    print(vertical_slice_umf_str)
    
    print(f"\nIn this formulation:")
    print(f" - The primary input variable is 'x'.")
    print(f" - The secondary input variable is 'u'.")
    print(f" - The mean of the primary UMF is c_upper = {c_upper}.")
    print(f" - The standard deviation of the primary UMF is σ_upper = {sigma_upper}.")
    print(f" - The spread of the vertical slice is controlled by β = {beta}.")

# Execute the function to print the formulation
generate_it3_formulation()