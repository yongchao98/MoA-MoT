def solve_superconductor_magnetization():
    """
    This script derives and prints the analytical expression for the initial
    magnetization curve of a superconducting bar in a transverse field.
    """
    print("Derivation of the Initial Magnetization Curve for a Superconducting Bar")
    print("=======================================================================")
    print("We consider a long superconducting bar with a rectangular cross-section (-a <= x <= a)")
    print("in an applied magnetic field H along the y-axis. The critical-state model is used,")
    print("assuming a constant critical current density Jc and the limit b << a.\n")

    print("Step 1: Relate Penetration Depth (d) to Applied Field (H)")
    print("-----------------------------------------------------------")
    print("From Ampere's law (dHy/dx = Jz) and boundary conditions (Hy(a)=H, Hy(a-d)=0),")
    print("we integrate to find the relationship between the field H and the penetration depth d.")
    print("The result is: d = H / Jc")
    print("This is valid for H < Hp, where Hp = a*Jc is the full penetration field.\n")

    print("Step 2: Calculate the Magnetic Moment (m_y)")
    print("--------------------------------------------")
    print("The magnetic moment per unit length (m_y) is found by integrating the contribution")
    print("from the shielding currents Jz over the cross-section A: m_y = integral(-x * Jz) dA.")
    print("Evaluating this integral gives:")
    print("m_y = -2b * Jc * d * (2a - d)\n")

    print("Step 3: Derive the Magnetization (M)")
    print("--------------------------------------")
    print("Magnetization M is the moment per unit volume (V = 4ab per unit length).")
    print("M = m_y / V = [-2b * Jc * d * (2a - d)] / (4ab)")
    print("This simplifies to: M = (-Jc / (2a)) * d * (2a - d)\n")

    print("Step 4: Final Expression for M(H)")
    print("-----------------------------------")
    print("Finally, we substitute d = H / Jc into the expression for M:")
    print("M = (-Jc / (2a)) * (H / Jc) * (2a - H / Jc)")
    print("Simplifying this gives the final analytical expression for the initial magnetization curve:\n")

    # Define variable names for the final equation printout
    M_var = "M"
    H_var = "H"
    a_var = "a"
    Jc_var = "Jc"
    one = "1"
    two = "2"

    # Print the final equation in a formatted box
    print("*******************************************************************")
    print("*                                                                 *")
    print(f"*    {M_var}({H_var}) = -{H_var} * ( {one} - {H_var} / ( {two} * {a_var} * {Jc_var} ) )                *")
    print("*                                                                 *")
    print(f"*    Valid for 0 <= {H_var} <= {H_var}p, where {H_var}p = {a_var} * {Jc_var}             *")
    print("*                                                                 *")
    print("*******************************************************************\n")

    print("Where:")
    print(f"  {M_var} = Magnetization (in A/m)")
    print(f"  {H_var} = Applied magnetic field (in A/m)")
    print(f"  {a_var} = Half-width of the bar (in m)")
    print(f"  {Jc_var} = Critical current density (in A/m^2)")

solve_superconductor_magnetization()
<<<M(H) = -H * (1 - H / (2 * a * Jc))>>>