import math

def display_force_law():
    """
    This function derives and displays the force law for a thermally isolated
    freely jointed chain polymer.
    """

    # The force law F(x) is the internal force of attraction between the polymer ends.
    # Based on the derivation in the microcanonical ensemble, the force is given by
    # the following equation. We will construct this equation string for display.

    # The numerical constants in the equation are:
    coefficient_const = 2
    power_const = 2

    # The equation shows the force F as a function of extension x.
    # It depends on the initial kinetic energy E(0), the number of segments n,
    # and the length of each segment l.
    
    force_law_expression = (
        f"F(x) = - ( ({coefficient_const} * E(0) * x) / (n**{power_const} * l**{power_const}) ) * "
        f"exp( x**{power_const} / (n**{power_const} * l**{power_const}) )"
    )

    print("The force law for the thermally isolated polymer is:")
    print(force_law_expression)
    print("\nWhere:")
    print("  F(x)  = The force of attraction between the polymer ends.")
    print("  x     = The separation distance between the polymer ends.")
    print("  E(0)  = The kinetic energy of the polymer at zero extension (x=0).")
    print("  n     = The number of segments in the polymer chain (assumed to be large).")
    print("  l     = The length of each segment.")
    print("  exp() = The exponential function, e^().")

# Execute the function to display the result.
display_force_law()