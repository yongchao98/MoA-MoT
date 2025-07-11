def print_force_law_isolated_polymer():
    """
    This function prints the derived force law for a thermally isolated
    freely jointed chain polymer, clearly showing all its components.
    """
    # Define the symbols used in the formula for clarity.
    F_x = "F(x)"
    x = "x"
    l = "l"
    n = "n"
    E0 = "E(0)"

    # Construct the final equation string.
    # The force is attractive (negative sign).
    # The coefficient part of the formula.
    coefficient = f"(2 * {E0} * {x} / ({n}**2 * {l}**2))"
    # The exponential part of the formula.
    exponent = f"({x}**2 / ({n}**2 * {l}**2))"
    exponential_term = f"exp({exponent})"
    
    final_equation = f"{F_x} = - {coefficient} * {exponential_term}"

    print("The force law for a thermally isolated freely jointed chain is:")
    print(final_equation)
    print("\nWhere:")
    print(f"  {F_x}: Force between the polymer ends at extension x.")
    print(f"  {x}: Separation of the polymer ends.")
    print(f"  {l}: Length of a single strut in the chain (ell).")
    print(f"  {n}: Number of mass points in the chain.")
    print(f"  {E0}: Kinetic energy of the polymer at zero extension (x=0).")
    print("  exp(): The exponential function e^().")

print_force_law_isolated_polymer()