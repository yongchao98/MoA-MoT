def generate_mori_tanaka_expression():
    """
    This function generates and prints the expression for the effective
    elastic moduli (C) in the Mori-Tanaka model.
    """

    # Define the symbols for the quantities as strings
    C = "C"
    I = "I"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"

    # Assemble the expression for C.
    # Note:
    # '+' and '-' represent tensor addition/subtraction.
    # '*' represents scalar multiplication.
    # ':' represents the double-dot product (tensor contraction on two indices).
    # '(...)^-1' represents the inverse of a tensor.
    
    print("The expression for the effective average elastic moduli C is:")
    
    # Printing each symbol in the final equation as requested
    # The f-string constructs the complete equation from its constituent parts.
    print(f"{C} = {Cm} + {Vf} * ({Cf} - {Cm}) : {A} : ({Vm} * {I} + {Vf} * {A})^-1")

# Execute the function
generate_mori_tanaka_expression()