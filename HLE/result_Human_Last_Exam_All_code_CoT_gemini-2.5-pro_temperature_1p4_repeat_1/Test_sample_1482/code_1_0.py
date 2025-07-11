import math

def generate_it3_vertical_slice_formula():
    """
    Generates and prints the mathematical formulation for the vertical cross-section
    of a Gaussian-based Interval Type-3 Membership Function (IT3 MF).
    """
    # Define the parameters for a concrete example of a Gaussian IT3 MF
    c = 5.0              # Center of the primary Gaussian
    underline_sigma = 1.0  # Lower bound of the primary standard deviation
    overline_sigma = 2.0   # Upper bound of the primary standard deviation
    k_u = 0.3            # Proportionality constant for the vertical slice's spread

    # Pre-calculate denominators for clarity in the final formula
    # Denominator term is 2 * sigma^2
    denom_upper = 2.0 * (overline_sigma**2)
    denom_lower = 2.0 * (underline_sigma**2)

    # The formula represents the upper bound of the membership grade for the
    # vertical slice at a specific primary input 'x'. This is a function of 'x'
    # and the secondary variable 'u'.
    #
    # The structure is: overline_mu(x, u) = exp(-0.5 * ((u - mean) / std_dev)**2)
    #
    # 1. The 'mean' is the upper bound of the primary IT2 MF:
    #    mean = exp(-((x - c)**2) / (2 * overline_sigma**2))
    #
    # 2. The 'std_dev' is proportional to the primary IT2 MF's uncertainty region:
    #    std_dev = k_u * (upper_bound_primary_MF - lower_bound_primary_MF)
    #

    # Construct the final formula string with all numbers explicitly included.
    formula = (
        f"overline_mu(x, u) = exp(-0.5 * ((u - exp(-(x - {c})**2 / {denom_upper})) / "
        f"({k_u} * (exp(-(x - {c})**2 / {denom_upper}) - exp(-(x - {c})**2 / {denom_lower}))))**2)"
    )
    
    print("The mathematical formulation for the upper bound of the vertical cross-section of a Gaussian IT3 MF is:")
    print(formula)

# Execute the function to print the formula
generate_it3_vertical_slice_formula()