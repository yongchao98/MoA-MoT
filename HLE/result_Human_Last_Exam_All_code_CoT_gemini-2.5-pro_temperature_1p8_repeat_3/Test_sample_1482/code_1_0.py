import math

def calculate_it3_mf_vertical_slice():
    """
    Calculates and displays the formulation for the upper bound of a vertical
    cross-section of an Interval Type-3 Membership Function (IT3 MF) based on
    Gaussian paradigms.
    """
    # --- Parameters ---
    # Primary input variable
    x = 4.5
    # Secondary input variable (determines the level of the slice)
    u = 0.6
    # Center of the Gaussian functions
    c = 5.0
    # Standard deviation for the "upper-upper" bounding Gaussian function.
    # This corresponds to the wider Gaussian.
    sigma_upper = 2.0
    # Standard deviation for the "lower-upper" bounding Gaussian function.
    # This corresponds to the narrower Gaussian.
    sigma_lower = 1.0

    # Ensure the constraint sigma_upper >= sigma_lower is met for a valid FOU
    if sigma_upper < sigma_lower:
        print("Error: sigma_upper must be greater than or equal to sigma_lower.")
        return

    # --- Formulation ---
    # The general formula for the upper bound of the vertical cross-section is:
    # ū_A(x, u) = u * ū_ūA(x) + (1 - u) * u̲_ūA(x)
    # where ū_ūA(x) and u̲_ūA(x) are the bounding Gaussian functions.
    # ū_ūA(x) = exp(-0.5 * ((x - c) / σ_upper)^2)
    # u̲_ūA(x) = exp(-0.5 * ((x - c) / σ_lower)^2)

    # --- Calculation ---
    # Calculate the value of the "upper-upper" bounding MF
    mf_upper_bound = math.exp(-0.5 * ((x - c) / sigma_upper)**2)

    # Calculate the value of the "lower-upper" bounding MF
    mf_lower_bound = math.exp(-0.5 * ((x - c) / sigma_lower)**2)

    # Calculate the final value of the vertical slice's upper bound
    result = u * mf_upper_bound + (1 - u) * mf_lower_bound

    # --- Output ---
    print("This script calculates the upper bound of a vertical cross-section for an IT3 MF.")
    print("--------------------------------------------------------------------------------\n")
    print("Mathematical Formulation:")
    print("overline_mu_A(x, u) = u * exp(-0.5 * ((x - c) / sigma_upper)^2) + (1 - u) * exp(-0.5 * ((x - c) / sigma_lower)^2)\n")
    print("Calculation with specific values:")

    # Printing the equation with all the numbers
    equation_str = (
        f"overline_mu_A({x}, {u}) = {u} * exp(-0.5 * (({x} - {c}) / {sigma_upper})^2) "
        f"+ (1 - {u}) * exp(-0.5 * (({x} - {c}) / {sigma_lower})^2)"
    )
    print(equation_str)

    # Show intermediate calculations for clarity
    calculation_str = (
        f"                   = {u} * {mf_upper_bound:.4f} + {1-u:.1f} * {mf_lower_bound:.4f}"
    )
    print(calculation_str)
    
    print(f"\nFinal Result: {result}")


if __name__ == "__main__":
    calculate_it3_mf_vertical_slice()