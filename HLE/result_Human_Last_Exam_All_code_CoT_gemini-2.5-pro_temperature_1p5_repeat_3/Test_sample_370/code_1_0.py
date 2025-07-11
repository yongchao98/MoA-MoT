def calculate_cross_section():
    """
    This function prints the step-by-step derivation of the total cross section
    for fermion-fermion scattering mediated by a pseudoscalar boson in the high-energy limit.
    """

    print("Calculation of the total cross section for psi + psi -> psi + psi\n")

    # Step 1: Explain the general formula and assumptions
    print("Step 1: General Formulas and Assumptions")
    print("-----------------------------------------")
    print("The total cross section sigma is found by integrating the differential cross section:")
    print("sigma = (1/2) * Integral( d(sigma)/d(Omega) d(Omega) )")
    print("The factor 1/2 accounts for the two identical fermions in the final state.")
    print("The differential cross section in the center-of-mass (CM) frame is:")
    print("d(sigma)/d(Omega) = (1 / (64 * pi**2 * s)) * overline|M|^2")
    print("Here, 's' is the Mandelstam variable s = (2*E)**2 = 4*E**2, and overline|M|^2 is the spin-averaged matrix element squared.")
    print("We assume the high-energy limit: fermion mass m=0 and CM energy E >> scalar mass M.\n")

    # Step 2: The matrix element
    print("Step 2: Spin-Averaged Matrix Element Squared")
    print("---------------------------------------------")
    print("The matrix element includes t-channel and u-channel diagrams. In the m=0 limit, the spin-averaged matrix element squared is:")
    print("overline|M|^2 = g**4 * ( t**2/(t - M**2)**2 + u**2/(u - M**2)**2 + (t*u)/((t - M**2)*(u - M**2)) )")
    print("where t and u are the other Mandelstam variables.\n")

    # Step 3: High-energy limit simplification
    print("Step 3: High-Energy Limit Simplification and Integration")
    print("---------------------------------------------------------")
    print("In the high-energy limit (E >> M), for most scattering angles, |t| >> M**2 and |u| >> M**2.")
    print("The expression for overline|M|^2 simplifies to a constant value:")
    print("overline|M|^2 -> g**4 * (t**2/t**2 + u**2/u**2 + (t*u)/(t*u)) = 3 * g**4")
    print("\nWe can now perform the integration:")
    print("sigma = (1/2) * Integral_d(phi) * Integral_d(cos(theta)) [ (1 / (64 * pi**2 * (4*E**2))) * (3 * g**4) ]")
    print("The integral over the solid angle gives 4*pi.")
    print("sigma = (1/2) * (4 * pi) * (1 / (256 * pi**2 * E**2)) * (3 * g**4)")
    print("sigma = (2 * pi / (256 * pi**2 * E**2)) * (3 * g**4)")
    print("sigma = (6 * pi * g**4) / (256 * pi**2 * E**2)\n")

    # Step 4: The final result
    print("Step 4: Final Result for the Total Cross Section")
    print("--------------------------------------------------")
    print("Simplifying the expression gives the final result for the total cross section:")
    print("Initial expression: sigma = (6 * pi * g**4) / (256 * pi**2 * E**2)")
    print("Final simplified equation: sigma = (3 * g**4) / (128 * pi * E**2)")

    # Print the equation with numbers highlighted as requested
    print("\nFinal equation with numerical factors shown explicitly:")
    print("sigma = (3 * g**4) / (128 * pi * E**2)")

calculate_cross_section()