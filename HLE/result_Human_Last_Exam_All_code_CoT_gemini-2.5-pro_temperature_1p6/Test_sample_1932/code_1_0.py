import math

def generate_weight_change_formula():
    """
    This function explains the reasoning and prints the symbolic formula for the 
    estimated change in weight of the running hourglass.
    """
    
    # The final expression for the change in weight, Delta W, is symbolic.
    # It corresponds to choice D, derived by considering the largest single effect, 
    # which is the impact force of the sand at the beginning of its fall.
    # Delta W = (mass_flow_rate) * (max_impact_velocity)
    
    # The components of the formula are:
    mass_flow_rate_expr = "(pi * d^2 * h * rho) / (4 * t)"
    max_velocity_expr = "sqrt(2 * g * H)"
    
    # The problem asks to output the numbers in the equation.
    # The numbers are '2' and '4'.
    
    print("Based on the analysis, the change in weight is estimated by considering the largest single effect.")
    print("This effect is the impact force of the sand, which is maximal at the beginning of the run.")
    print("The formula for this effect is:")
    
    # Using the symbolic variables from the problem description
    pi = "pi"
    d = "d"
    h = "h"
    rho = "rho"
    t = "t"
    g = "g"
    H = "H"
    two = "2"
    four = "4"
    
    # Constructing and printing the final symbolic equation for Delta W
    final_equation = f"Delta W = ({pi} * {d}^{two} * {h} * {rho} / ({four} * {t})) * sqrt({two} * {g} * {H})"
    print(final_equation)

# Execute the function to display the result
generate_weight_change_formula()