def solve_magnetization():
    """
    This function derives and prints the analytical expression for the initial
    magnetization curve of a superconducting bar based on the Bean critical-state model.
    """

    print("Step-by-step derivation of the initial magnetization curve M(H):")
    print("-" * 60)

    # Step 1: Explain the model
    print("Step 1: Define the Physical Model")
    print("We model the wide, thin bar (b << a) as an infinite superconducting slab of thickness 2a.")
    print("The magnetic field H is applied parallel to the slab's surfaces and penetrates from the sides at x = +/-a.")
    print("We use the Bean critical-state model, where the induced shielding current has a constant density, Jc.")
    print("\n")

    # Step 2: Field Penetration
    print("Step 2: Field Penetration in the Critical-State Model")
    print("As the external field H increases from zero, shielding currents penetrate from the surfaces to a depth 'd'.")
    print("The penetration depth 'd' is linearly related to the applied field H:")
    print("H = Jc * d")
    print("This is valid until the slab is fully penetrated (d=a).")
    print("\n")

    # Step 3: Magnetic Moment and Magnetization
    print("Step 3: Calculate the Magnetization M")
    print("The shielding currents create a magnetic moment opposing the applied field. The resulting magnetization M can be found by integrating the current distribution.")
    print("For a given penetration depth d, the magnetization M is:")
    print("M = -(Jc / (2*a)) * (2*a*d - d^2)")
    print("\n")

    # Step 4: Express M in terms of H
    print("Step 4: Express Magnetization as a Function of H")
    print("Substitute the relation d = H / Jc from Step 2 into the expression for M:")
    print("M = -(Jc / (2*a)) * (2*a*(H/Jc) - (H/Jc)^2)")
    print("After simplifying, we get:")
    print("M = -(H - H^2 / (2 * a * Jc))")
    print("Factoring out -H gives the final form of the equation.")
    print("\n")

    # Final Answer
    print("Final Analytical Expression:")
    print("-" * 60)
    print("The analytical expression for the initial magnetization curve is:")

    a = 'a'
    Jc = 'Jc'
    H = 'H'

    # The numbers in the equation are 1 and 2.
    final_equation_str = f"M({H}) = -{H} * ({1} - {H} / ({2}*{a}*{Jc}))"
    print(final_equation_str)

    print("\nThis expression is valid for an increasing field H from 0 up to the full penetration field, Hp.")
    print("The full penetration field Hp is the field at which d=a:")
    print(f"Hp = {Jc} * {a}")
    print("In terms of Hp, the expression is:")
    Hp = "Hp"
    final_equation_Hp_str = f"M({H}) = -{H} * ({1} - {H} / ({2}*{Hp}))   for 0 <= H <= {Hp}"
    print(final_equation_Hp_str)

solve_magnetization()
<<<M(H) = -H * (1 - H / (2*a*Jc))>>>