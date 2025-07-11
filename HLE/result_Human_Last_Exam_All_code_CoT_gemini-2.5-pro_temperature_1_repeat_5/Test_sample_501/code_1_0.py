import math

def print_force_law():
    """
    This function prints the derived force law for a thermally isolated polymer chain.
    The force F(x) is given as a function of the extension x, the number of links n,
    the length of each link l, and the initial kinetic energy at zero extension E(0).
    """

    # The equation has the structure: F(x) = A * x * exp(B * x^2)
    # We will define the components A and B based on the derivation.

    # Coefficient for the linear term
    coeff_A_numerator = 3
    coeff_A_denominator = "n^2 * l^2"

    # Coefficient for the term in the exponent
    coeff_B_numerator = -3
    coeff_B_denominator = "2 * n^2 * l^2"

    print("The derived force law, F(x), for the thermally isolated polymer is:")
    print("F(x) = (C1 * E(0) * x / (n^2 * l^2)) * exp(C2 * x^2 / (n^2 * l^2))")
    print("\nWhere E(0) is the kinetic energy at zero extension (x=0).")
    print("The numerical constants in the equation are:")

    print(f"C1 = {coeff_A_numerator}")
    c2_ratio = coeff_B_numerator / 2
    print(f"C2 = {coeff_B_numerator} / 2 = {c2_ratio}")

    # Reconstruct and print the final, complete formula string
    final_formula = f"F(x) = ({coeff_A_numerator} * E(0) * x / ({coeff_A_denominator})) * exp({coeff_B_numerator} * x^2 / ({coeff_B_denominator}))"
    print("\nTherefore, the final equation is:")
    print(final_formula)

if __name__ == "__main__":
    print_force_law()