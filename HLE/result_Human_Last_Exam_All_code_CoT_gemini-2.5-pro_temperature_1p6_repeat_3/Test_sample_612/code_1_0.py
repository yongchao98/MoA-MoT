def get_mori_tanaka_expression():
    """
    This function constructs and prints the symbolic expression for the
    effective average elastic moduli (C) using the Mori-Tanaka model.
    """

    # Define the symbols for the quantities as strings
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    I = "I"
    A = "A"

    # Construct the components of the formula as strings
    # In tensor algebra, (X)^(-1) is often written as inv(X)
    stiffness_difference = f"({Cf} - {Cm})"
    strain_relation_term = f"((1 - {Vf})*{I} + {Vf}*{A})"
    inverse_term = f"inv({strain_relation_term})"

    # Assemble the final expression for C
    # The formula is C = Cm + Vf * (Cf - Cm) : A : [ ((1 - Vf)*I + Vf*A) ]^-1
    # We use ":" to represent the double-dot product.
    final_expression = f"{Cm} + {Vf}*{stiffness_difference} : {A} : {inverse_term}"

    print("The expression for the effective average elastic moduli tensor C is:")
    print(f"\nC = {final_expression}\n")
    print("Where each symbol represents:")
    print(f" C: The effective average elastic moduli tensor")
    print(f" {Cm}: The fourth-order elasticity tensor of the polymer matrix")
    print(f" {Cf}: The fourth-order elasticity tensor of the fiber")
    print(f" {Vf}: The volume fraction of the fibers (so Vm = 1 - Vf)")
    print(f" {I}: The fourth-order identity tensor")
    print(f" {A}: The Eshelby strain-concentration tensor")
    print("\nAnd the operators are:")
    print(" :   Represents the double-dot product between tensors.")
    print(" *   Represents scalar multiplication.")
    print("inv(...) Represents the inverse of the enclosed fourth-order tensor.")

# Execute the function to print the result
get_mori_tanaka_expression()