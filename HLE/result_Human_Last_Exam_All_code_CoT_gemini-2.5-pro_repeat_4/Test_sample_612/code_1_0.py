def generate_mori_tanaka_expression():
    """
    This function generates and prints the expression for the effective average
    elastic moduli (C) in the Mori-Tanaka model for composite materials.
    """
    # Define the symbols for the quantities involved
    I = "I"      # The fourth-order identity tensor
    Cf = "Cf"    # The fourth-order elasticity tensor of the fiber
    Cm = "Cm"    # The fourth-order elasticity tensor of the polymer matrix
    Vf = "Vf"    # The volume fraction of the fibers
    Vm = "Vm"    # The volume fraction of the matrix
    A = "A"      # The Eshelby strain-concentration tensor
    C = "C"      # The effective average elastic moduli tensor

    # Construct the two main parts of the expression
    # Note: '*' denotes tensor multiplication and '(...)^-1' denotes inversion.
    
    # First part: Volume average of the 'effective' phase stiffnesses
    # This is (Vf * Cf * A + Vm * Cm)
    term1 = f"({Vf} * {Cf} * {A} + {Vm} * {Cm})"
    
    # Second part: The inverse of the volume-averaged strain concentration
    # This is (Vf * A + Vm * I)^-1
    term2 = f"({Vf} * {A} + {Vm} * {I})**-1"

    # Combine the parts to form the final expression for C
    expression = f"{C} = {term1} * {term2}"

    # Print the final expression, showing each symbol in the equation
    print("The Mori-Tanaka expression for the effective average elastic moduli C is:")
    print(expression)

if __name__ == "__main__":
    generate_mori_tanaka_expression()