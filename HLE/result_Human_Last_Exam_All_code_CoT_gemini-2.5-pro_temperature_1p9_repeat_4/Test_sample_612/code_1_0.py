def get_mori_tanaka_expression():
    """
    This function prints the symbolic expression for the effective elastic moduli (C)
    of a fiber-reinforced composite using the Mori-Tanaka model.
    It also lists the symbols used in the equation as requested.
    """

    # Define the symbolic variables to build the final expression string.
    C = "C"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"
    I = "I"

    # Construct the final equation string. In this text representation,
    # '*' denotes tensor multiplication/contraction and '**-1' denotes inversion.
    final_equation = f"{C} = {Cm} + {Vf} * ({Cf} - {Cm}) * {A} * ({Vm} * {I} + {Vf} * {A})**-1"

    print("The expression for the effective average elastic moduli C is:")
    print(final_equation)

    print("\nThe individual symbols in the final equation are:")
    print(C)
    print(Cm)
    print(Vf)
    print(Cf)
    print(A)
    print(Vm)
    print(I)

# Execute the function to display the answer.
get_mori_tanaka_expression()