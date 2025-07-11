def display_mori_tanaka_expression():
    """
    This function prints the expression for the effective average elastic moduli (C)
    in the Mori-Tanaka model using the provided symbolic quantities.
    """
    # Define the symbols used in the equation as string variables.
    symbol_C = "C"
    symbol_I = "I"
    symbol_Cf = "Cf"
    symbol_Cm = "Cm"
    symbol_Vf = "Vf"
    symbol_A = "A"

    print("The expression for the effective average elastic moduli (C) in the Mori–Tanaka model is given by the following equation:")
    print("-" * 80)

    # We use 'end=""' to print all parts on the same line, which constructs the equation step-by-step.
    # This fulfills the request to output each symbol in the final equation.
    print(f"{symbol_C} = ", end="")
    print(f"{symbol_Cm} + ", end="")
    print(f"{symbol_Vf} * ", end="")
    print(f"({symbol_Cf} - {symbol_Cm}) * ", end="")
    print(f"{symbol_A} * ", end="")
    print(f"[(1 - {symbol_Vf}) * {symbol_I} + {symbol_Vf} * {symbol_A}]^-1")
    print("-" * 80)
    print("Where:")
    print(f"  C: The effective average elastic moduli tensor of the composite")
    print(f"  Cm: The elasticity tensor of the matrix")
    print(f"  Cf: The elasticity tensor of the fiber")
    print(f"  Vf: The volume fraction of the fiber")
    print(f"  I:  The fourth-order identity tensor")
    print(f"  A:  The Eshelby strain-concentration tensor")
    print(f"  (...)^-1: Denotes the inverse of the tensor")

# Execute the function to print the expression.
display_mori_tanaka_expression()
