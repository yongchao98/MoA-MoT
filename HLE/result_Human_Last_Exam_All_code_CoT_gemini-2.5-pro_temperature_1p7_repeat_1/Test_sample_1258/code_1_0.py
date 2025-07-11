import math

def print_demagnetizing_factor_expression():
    """
    Prints the analytical expression for the fluxmetric demagnetizing factor
    for a magnetic cylinder.
    """
    
    # The expression is broken down into parts for clarity.
    
    # Header
    print("The analytical expression for the fluxmetric demagnetizing factor (N_f) is:")
    print("-" * 70)
    
    # Main formula structure
    print("N_f = (4 / (3 * pi * g**2)) * [Term_A - Term_B - Term_C]")
    print("") # Blank line for spacing
    
    # Definitions of terms
    print("Where:")
    print("  Term_A = sqrt(g**2 + 1) * (g**2 + 2) * E(k)")
    print("  Term_B = (1 / sqrt(g**2 + 1)) * F(k)")
    print("  Term_C = g**2 + 1/2")
    print("")
    
    # Full expression in one line
    print("Putting it all together:")
    print("N_f = (4 / (3 * pi * g**2)) * (sqrt(g**2 + 1)*(g**2 + 2)*E(k) - (1/sqrt(g**2 + 1))*F(k) - g**2 - 0.5)")
    print("-" * 70)
    
    # Definitions of variables
    print("And the variables are defined as:")
    print("  g = Length-to-Diameter ratio of the cylinder (L/D)")
    print("  k = Modulus of the elliptic integrals, defined by k**2 = 1 / (1 + g**2 / 4)")
    print("  F(k) = Complete elliptic integral of the first kind with modulus k")
    print("  E(k) = Complete elliptic integral of the second kind with modulus k")
    print("  pi = The mathematical constant " + str(math.pi))
    print("  The numbers in the equation are: 4, 3, 1, 2, 0.5")

if __name__ == '__main__':
    print_demagnetizing_factor_expression()
