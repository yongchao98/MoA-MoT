def get_mori_tanaka_expression():
    """
    This function constructs and prints the symbolic expression for the
    effective average elastic moduli (C) using the Mori-Tanaka model.
    """
    # Define the symbols as strings for the equation
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vm = "Vm"
    Vf = "Vf"
    A = "A"
    I = "I"

    # Construct the terms of the equation as strings
    # The final expression is: C = Cm + Vf * (Cf - Cm) * A * inv(Vm*I + Vf*A)
    # We use '**-1' to denote the tensor inverse.
    # Tensor multiplication is represented by '*'.

    term1 = Cm
    term2 = f"{Vf} * ({Cf} - {Cm})"
    term3 = A
    inverse_term = f"({Vm}*{I} + {Vf}*{A})**-1"

    # Assemble the final equation string
    final_equation = f"{C} = {term1} + {term2} * {term3} * {inverse_term}"

    print("The expression for the effective average elastic moduli C in the Mori-Tanaka model is:")
    print(final_equation)
    print("\nWhere:")
    print(f"  {C}: Effective average elastic moduli of the composite")
    print(f"  {Cm}: Elasticity tensor of the matrix")
    print(f"  {Cf}: Elasticity tensor of the fiber")
    print(f"  {Vm}: Volume fraction of the matrix")
    print(f"  {Vf}: Volume fraction of the fiber")
    print(f"  {A}: Strain-concentration tensor")
    print(f"  {I}: Fourth-order identity tensor")
    print("  '*' denotes tensor multiplication and '**-1' denotes tensor inversion.")


if __name__ == '__main__':
    get_mori_tanaka_expression()