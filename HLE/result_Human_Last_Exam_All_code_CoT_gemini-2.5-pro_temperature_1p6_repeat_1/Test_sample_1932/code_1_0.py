import math

def display_weight_change_formula():
    """
    This function prints the formula for the estimated change in weight of a running hourglass.
    The chosen formula corresponds to option A from the problem description.
    """
    # The formula is Delta_W = (mass_flow_rate) * (impact_velocity)
    # mass_flow_rate = (pi * d**2 * h * rho) / (4 * t)
    # impact_velocity = sqrt(2 * g * (H - h))

    # The formula contains the numbers 4 and 2, which will be part of the printed output.
    numerator = "pi * d**2 * h * rho"
    denominator = "4 * t"
    sqrt_content = "2 * g * (H - h)"
    
    # We construct and print the formula as a readable string.
    print("The expression for the estimated change in weight (Delta W) is:")
    print(f"Delta W = ({numerator} / {denominator}) * sqrt({sqrt_content})")

display_weight_change_formula()