import math

def generate_it3_formulation():
    """
    Generates and prints the mathematical formulation for the vertical
    cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).

    This formulation describes the uncertainty bounds for a fixed primary input 'x'
    and a secondary input 'u'.
    """

    # --- Parameters for the IT3 Gaussian Membership Function ---
    # These values are chosen for illustrative purposes.

    # Center (mean) of the Gaussian function for the primary variable 'x'.
    center_x = 5.0

    # Standard deviation for the upper bound of the primary membership function of 'x'.
    # This determines the 'width' of the main bell curve.
    sigma_x_upper = 2.0

    # Standard deviation for the Gaussian in the secondary 'u' dimension.
    # This models the 'blurriness' or uncertainty in the vertical slice.
    sigma_u_upper = 0.3

    # --- Construct the Mathematical Formulation ---
    # The formulation represents the Upper Membership Function (UMF) of the
    # Type-2 Fuzzy Set that constitutes the vertical cross-section of the IT3 MF.
    # It is denoted as μ_upper(x, u).
    # The structure is a Gaussian in 'u', where its mean is another Gaussian in 'x'.

    # The equation string uses f-string formatting to insert the parameters.
    # The format exp(-((z - mean)/std_dev)**2) is a common convention for Gaussian MFs.
    equation_str = (
        f"μ_upper(x, u) = exp(-((u - exp(-((x - {center_x}) / {sigma_x_upper})**2)) / {sigma_u_upper})**2)"
    )

    # --- Print the Final Output ---
    print("The mathematical formulation for the upper bound of a vertical cross-section of a Gaussian IT3 MF is:")
    print(equation_str)

if __name__ == "__main__":
    generate_it3_formulation()