def print_mori_tanaka_equation():
    """
    This function prints the expression for the effective average elastic moduli (C)
    in the Mori-Tanaka model for a fiber-reinforced composite.
    """
    # Define the symbols as string variables for clarity.
    # These represent fourth-order tensors and scalar volume fractions.
    C = "C"    # Effective average elastic moduli
    Cm = "Cm"  # Elasticity tensor of the matrix
    Cf = "Cf"  # Elasticity tensor of the fiber
    Vf = "Vf"  # Volume fraction of the fiber
    Vm = "Vm"  # Volume fraction of the matrix
    A = "A"    # Eshelby strain-concentration tensor
    I = "I"    # Fourth-order identity tensor

    # Construct the Mori-Tanaka equation as a formatted string.
    # The ':' operator is used here to represent the double-dot product in tensor notation.
    # The notation '[... ]^-1' represents the inverse of the enclosed fourth-order tensor.
    equation = f"{C} = {Cm} + {Vf} * ({Cf} - {Cm}) : {A} : [{Vm} * {I} + {Vf} * {A}]^-1"

    # Define the relationship between the volume fractions
    volume_relation = f"where {Vm} + {Vf} = 1."

    # Print the final expression and the associated relation
    print("The expression for the effective average elastic moduli (C) is:")
    print(equation)
    print(volume_relation)

# Execute the function to display the result
print_mori_tanaka_equation()