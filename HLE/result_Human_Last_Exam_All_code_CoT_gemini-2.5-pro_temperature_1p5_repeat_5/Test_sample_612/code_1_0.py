def get_mori_tanaka_expression():
    """
    This function prints the symbolic expression for the effective average elastic moduli (C)
    in a fiber-reinforced composite, according to the Mori-Tanaka model.
    """

    # Define the symbolic variables as strings for printing
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vf = "Vf"
    Vm = "Vm"
    I = "I"
    A = "A"

    # The Mori-Tanaka equation for the effective stiffness tensor C is constructed.
    # In tensor mechanics notation, multiplication '*' implies a tensor product (e.g., a double dot product),
    # and '^-1' denotes a tensor inversion.
    
    # Print each part of the final equation to fulfill the user request
    print("The Mori-Tanaka expression for the effective elastic moduli C is:")
    print(f"{C} = {Cm} + {Vf} * ({Cf} - {Cm}) * {A} * [{Vm}*{I} + {Vf}*{A}]**-1")

get_mori_tanaka_expression()