import math

def display_mutual_inductance_change():
    """
    This function prints the derived expression for the change in mutual
    inductance per unit length between the two circuits when the
    concentrator is added.

    The expression is: Delta_M = (mu_0 * h^2) / (2 * pi * R1^2)

    Where:
    - mu_0 is the permeability of free space.
    - h is the separation distance between the wires in each circuit.
    - R1 is the inner radius of the concentrator shell.
    - pi is the mathematical constant.
    """

    # Symbolic components of the final equation
    mu_0_symbol = "mu_0"
    h_symbol = "h"
    R1_symbol = "R1"
    pi_symbol = "pi"

    # Numeric components of the final equation
    h_power = 2
    R1_power = 2
    denominator_coefficient = 2

    # Print the explanation and the final formula as a string
    print("The derived expression for the change in mutual inductance per unit length (Delta_M) is:")
    
    final_expression = (f"({mu_0_symbol} * {h_symbol}**{h_power}) / "
                        f"({denominator_coefficient} * {pi_symbol} * {R1_symbol}**{R1_power})")
    
    print("\nDelta_M = " + final_expression)
    
    # Per the instruction to output each number in the final equation,
    # we'll print them explicitly here.
    print("\nBreaking down the expression:")
    print(f"Numerator term: {mu_0_symbol} * {h_symbol}^{h_power}")
    print(f"Denominator term: {denominator_coefficient} * {pi_symbol} * {R1_symbol}^{R1_power}")
    print("\nNumeric constants in the formula:")
    print(f"Power of h: {h_power}")
    print(f"Coefficient in denominator: {denominator_coefficient}")
    print(f"Power of R1: {R1_power}")


if __name__ == "__main__":
    display_mutual_inductance_change()
