def generate_it3_vertical_slice_formulation():
    """
    Generates and prints the mathematical formulation for the vertical cross-section
    of a Gaussian Interval Type-3 Membership Function (IT3 MF).
    """

    # --- Parameters for the IT3 MF and the Vertical Slice ---

    # Center of the primary Gaussian membership function
    c = 5.0
    # The fixed value of the primary input 'x' for the vertical cross-section
    x_prime = 6.0

    # For the linear functions f(u) = a*u + b, which define how the Gaussian's
    # parameters (scale 's' and standard deviation 'sigma') change with 'u'.
    
    # Parameters for the upper scaling function: s_upper(u) = a_s_upper * u + b_s_upper
    # Represents s_upper changing from 1.0 (at u=0) to 1.0 (at u=1)
    a_s_upper = 0.0
    b_s_upper = 1.0

    # Parameters for the lower scaling function: s_lower(u) = a_s_lower * u + b_s_lower
    # Represents s_lower changing from 0.8 (at u=0) to 0.9 (at u=1)
    a_s_lower = 0.1
    b_s_lower = 0.8

    # Parameters for the upper std dev function: σ_upper(u) = a_sigma_upper * u + b_sigma_upper
    # Represents σ_upper changing from 2.0 (at u=0) to 2.5 (at u=1)
    a_sigma_upper = 0.5
    b_sigma_upper = 2.0

    # Parameters for the lower std dev function: σ_lower(u) = a_sigma_lower * u + b_sigma_lower
    # Represents σ_lower changing from 1.0 (at u=0) to 1.5 (at u=1)
    a_sigma_lower = 0.5
    b_sigma_lower = 1.0

    # --- Construct and Print the Formulation ---

    # Double curly braces {{ }} are used within the f-string to output a literal curly brace.
    upper_bound_equation = (f"μ_upper(u | x={x_prime}) = "
                            f"({a_s_upper} * u + {b_s_upper}) * "
                            f"exp(-0.5 * (({x_prime} - {c}) / ({a_sigma_upper} * u + {b_sigma_upper}))^2)")

    lower_bound_equation = (f"μ_lower(u | x={x_prime}) = "
                            f"({a_s_lower} * u + {b_s_lower}) * "
                            f"exp(-0.5 * (({x_prime} - {c}) / ({a_sigma_lower} * u + {b_sigma_lower}))^2)")

    print("The mathematical formulation for the vertical cross-section of a Gaussian Interval Type-3 MF")
    print(f"at a fixed primary input x' = {x_prime} is given by the following upper and lower bounds, as functions of u:")
    print("-" * 80)
    
    print("\n[Upper Bound Formulation]")
    print(upper_bound_equation)
    
    print("\n[Lower Bound Formulation]")
    print(lower_bound_equation)
    
    print("-" * 80)
    print("Where:")
    print(" -> 'u' is the secondary input variable (secondary membership grade), with u ∈ [0, 1].")
    print(" -> 'exp()' is the exponential function.")
    print(f" -> The formula represents the uncertainty bounds for a vertical slice taken at x' = {x_prime} from an IT3 MF centered at c = {c}.")

# Execute the function to print the formulation
generate_it3_vertical_slice_formulation()