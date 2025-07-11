def generate_mori_tanaka_expression():
    """
    This function generates and prints the symbolic expression for the
    effective average elastic moduli (C) in the Mori-Tanaka model.
    """
    # Define the symbols for the physical quantities
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vf = "Vf"
    I = "I"
    A = "A"

    # Assemble the components of the Mori-Tanaka equation
    # The model is C = Cm + Vf * (Cf - Cm) * A_MT,
    # where the Mori-Tanaka strain concentration tensor A_MT is A * [Vm*I + Vf*A]⁻¹
    # We substitute Vm with (1 - Vf) for an expression purely in terms of fiber volume fraction.

    part1 = Cm
    part2 = Vf
    part3 = f"({Cf} - {Cm})"
    part4 = A
    inverse_part = f"[(1 - {Vf}) * {I} + {Vf} * {A}]"

    # In tensor notation, multiplication is represented by operators like '*' or ':'.
    # We will use '*' as a general multiplication symbol for clarity.
    # The '⁻¹' symbol indicates the tensor inverse.
    final_expression = f"{C} = {part1} + {part2} * {part3} * {part4} * {inverse_part}⁻¹"

    # Print the final formatted equation
    print(final_expression)

# Execute the function to get the expression
generate_mori_tanaka_expression()