def solve_mutual_inductance_change():
    """
    This function prints the derived symbolic expression for the change in mutual
    inductance per unit length between the two circuits when the concentrators are added.
    """

    # Define the symbols used in the expression as strings for clear printing.
    mu_0 = "mu_0"
    pi = "pi"
    h = "h"
    d = "d"
    R1 = "R1"
    R2 = "R2"
    
    # Construct the final expression as a formatted string.
    # The numbers 2 and 1 from the derivation are explicitly included in the output.
    expression = f"({mu_0} * {h}**2 / (2 * {pi} * {d}**2)) * (({R2} / {R1})**2 - 1)"
    
    print("The expression for the change in mutual inductance per unit length, Delta_M, is:")
    print(f"Delta_M = {expression}")

solve_mutual_inductance_change()