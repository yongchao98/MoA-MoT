def solve_mori_tanaka_expression():
    """
    This function prints the symbolic expression for the effective average elastic moduli (C)
    in a fiber-reinforced composite based on the Mori-Tanaka model.

    The expression uses the following symbolic quantities:
    C:  Effective average elastic moduli tensor
    I:  Fourth-order identity tensor
    Cf: Fourth-order elasticity tensor of the fiber
    Cm: Fourth-order elasticity tensor of the polymer matrix
    Vf: Volume fraction of the fibers
    Vm: Volume fraction of the matrix
    A:  Eshelby strain-concentration tensor
    """

    # The string below represents the derived formula for C.
    # The colon ':' represents the double-dot product (a tensor contraction over two indices).
    # The function 'inv()' represents the tensor inverse operation.
    equation = "C = Cm + Vf * (Cf - Cm) : A : inv(Vm * I + Vf * A)"

    print("The expression for the effective average elastic moduli C is:")
    print(equation)
    print("\nWhere:")
    print(" C  = Effective average elastic moduli tensor")
    print(" I  = Fourth-order identity tensor")
    print(" Cf = Elasticity tensor of the fiber")
    print(" Cm = Elasticity tensor of the matrix")
    print(" Vf = Volume fraction of the fiber")
    print(" Vm = Volume fraction of the matrix (equal to 1 - Vf)")
    print(" A  = Eshelby strain-concentration tensor")
    print(" :  = Double-dot product")
    print("inv() = Tensor inverse")

solve_mori_tanaka_expression()