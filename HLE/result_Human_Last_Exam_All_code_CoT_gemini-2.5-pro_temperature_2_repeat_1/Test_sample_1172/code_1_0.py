import math

def display_inductance_change_formula():
    """
    This function prints the derived formula for the change in mutual inductance
    per unit length between the two circuits when concentrators are added.
    """
    
    # The derived change in mutual inductance per unit length is:
    # Delta_M' = - (mu_0 * h^2) / (2 * pi * d^2)
    # where mu_0 is the permeability of free space (4 * pi * 10^-7 H/m)
    
    # Substituting the value of mu_0 and simplifying:
    # Delta_M' = - (4 * pi * 10^-7 * h^2) / (2 * pi * d^2)
    # Delta_M' = - 2 * 10^-7 * (h^2 / d^2)
    
    coefficient = -2e-7

    print("The expression for the change in mutual inductance per unit length (Delta_M') is:")
    print("==================================================================================")
    print(f"Delta_M' = C * (h^2 / d^2)")
    print("\nwhere:")
    print("  - h: separation distance between wires in a circuit (in meters).")
    print("  - d: separation distance between the two circuits (in meters), with d >> h.")
    print("  - C: a constant coefficient.")
    print("\nThe values of the components in the equation are:")
    # Per the instruction to output each number in the final equation
    print(f"  - The numerical coefficient C = {coefficient:.1e} H/m.")
    print("  - The exponent for h is 2.")
    print("  - The exponent for d is 2.")
    print("\nSo the full expression is:")
    print(f"Delta_M' = ({coefficient:.1e}) * (h**2 / d**2)  (in Henries/meter)")

if __name__ == "__main__":
    display_inductance_change_formula()