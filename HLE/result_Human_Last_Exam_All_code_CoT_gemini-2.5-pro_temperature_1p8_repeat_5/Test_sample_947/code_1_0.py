def solve_superconductor_field():
    """
    This function prints the derived expression for the magnetic field
    around a superconducting strip under the specified conditions.
    """

    # Define the variables for clarity in the final expression
    H_z = "H_z(x, z)"
    Ha = "Ha"
    H0 = "H0"
    Jc = "Jc"
    d = "d"
    w = "w"
    x = "x"
    z = "z"
    pi = "pi"

    # Step 1: State the base formula for the total magnetic field
    print("The total magnetic field H_z is the sum of the applied field Ha and the induced field H_z_ind.")
    print(f"{H_z} = {Ha} + H_z_ind(x, z)\n")

    # Step 2: Provide the expression for the induced field
    print("The exact induced field H_z_ind for a single strip is:")
    print(f"H_z_ind = ({Jc}*{d} / (2*{pi})) * [arctan(2*{x}*{z} / ({x}^2 - {z}^2 - a^2)) - arctan(2*{x}*{z} / ({x}^2 - {z}^2 - {w}^2))]\n")

    # Step 3: Explain and apply the approximation for |x| >> a
    print(f"Under the approximation |{x}| >> a, we expand the first arctan term with respect to a^2:")
    print(f"arctan(2*{x}*{z} / ({x}^2 - {z}^2 - a^2)) approx = arctan(2*{x}*{z} / ({x}^2 - {z}^2)) + a^2 * (2*{x}*{z}) / ({x}^2 + {z}^2)^2\n")

    # Step 4: Substitute the approximation back into the formula for H_z
    print("Substituting this approximation, the total field becomes:")
    term1 = f"({Jc}*{d} / (2*{pi})) * [arctan(2*{x}*{z} / ({x}^2 - {z}^2)) - arctan(2*{x}*{z} / ({x}^2 - {z}^2 - {w}^2))]"
    term2 = f"({Jc}*{d} / (2*{pi})) * a^2 * (2*{x}*{z}) / ({x}^2 + {z}^2)^2"
    print(f"{H_z} approx = {Ha} + {term1} + {term2}\n")

    # Step 5: Relate 'a' to Ha using the given H0
    print("The flux front position 'a' is related to the applied field Ha by a = w * exp(-Ha / H0), where H0 = Jc*d/pi.")
    print(f"So, a^2 = {w}^2 * exp(-2 * {Ha} / {H0})\n")

    # Step 6: Present the final expression
    print("Substituting a^2, the final expression for the magnetic field is:")
    final_term1 = f"({H0}/2) * [arctan(2*{x}*{z} / ({x}^2 - {z}^2)) - arctan(2*{x}*{z} / ({x}^2 - {z}^2 - {w}^2))]"
    final_term2 = f"{H0} * ({x}*{z}*{w}^2) / ({x}^2 + {z}^2)^2 * exp(-2*{Ha}/{H0})"
    final_expression = f"{H_z} = {Ha} + {final_term1} + {final_term2}"
    
    print(final_expression)
    
    # Final answer as per requested format
    # Because the request asks for an expression, the final answer is the string of the formula.
    print("\n<<<FINAL EXPRESSION>>>")
    print(f"<<<{final_expression}>>>")

solve_superconductor_field()