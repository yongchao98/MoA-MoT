import sympy

def get_force_law():
    """
    This function calculates and prints the force law for a thermally isolated,
    freely jointed polymer chain.

    The derivation uses the microcanonical ensemble, where the total energy E
    is constant. The force of attraction between the ends of the polymer is
    derived from the principles of statistical mechanics.
    """
    
    # Define the symbols used in the equation
    E0 = sympy.Symbol("E(0)") # Kinetic energy at zero extension
    n = sympy.Symbol("n")     # Number of struts
    l = sympy.Symbol("ell")   # Length of each strut
    x = sympy.Symbol("x")     # Separation of the polymer ends

    # The constant factor in the force law
    constant_factor = 3
    
    # The derived expression for the force of attraction
    # F_attr = (3 * E(0) * x) / (n * (n + 1) * l**2)
    # For large n, this is approximately (3 * E(0) * x) / (n**2 * l**2)

    # We will print the more accurate formula using (n+1)*n
    print("The derived force law for the force of attraction between the polymer ends is:")
    print("\nF_attr = ({c} * {E0} * {x}) / ({n} * ({n} + 1) * {l}**2)".format(
        c=constant_factor,
        E0=E0,
        x=x,
        n=n,
        l=l
    ))
    
    print("\nWhere:")
    print(f"  F_attr: The force of attraction between the polymer ends.")
    print(f"  {constant_factor}: A numerical constant from the derivation.")
    print(f"  {E0}: The total kinetic energy of the polymer at zero extension.")
    print(f"  {x}: The separation of the polymer ends (assumed to be small).")
    print(f"  {n}: The number of links in the polymer chain (assumed to be large).")
    print(f"  {l}: The length of a single link.")

if __name__ == '__main__':
    get_force_law()
