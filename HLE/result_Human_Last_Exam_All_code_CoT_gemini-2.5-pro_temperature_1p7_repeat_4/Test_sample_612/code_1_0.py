import sympy

def mori_tanaka_equation():
    """
    This function prints the symbolic expression for the effective average elastic moduli (C)
    of a fiber-reinforced composite based on the Mori-Tanaka model.
    """

    # Define the symbols for the quantities involved.
    # Using sympy to represent them, though we only use their string representation here.
    C = sympy.Symbol("C")
    Cm = sympy.Symbol("Cm")
    Cf = sympy.Symbol("Cf")
    Vm = sympy.Symbol("Vm")
    Vf = sympy.Symbol("Vf")
    A = sympy.Symbol("A")
    I = sympy.Symbol("I")

    # The Mori-Tanaka model is based on the following derivation:
    # 1. Average stress: <sigma> = Vm*<sigma_m> + Vf*<sigma_f>
    # 2. Average strain: <epsilon> = Vm*<epsilon_m> + Vf*<epsilon_f>
    # 3. Constitutive laws: <sigma_m> = Cm*<epsilon_m>, <sigma_f> = Cf*<epsilon_f>
    # 4. Mori-Tanaka hypothesis: <epsilon_f> = A*<epsilon_m>
    # 5. From (2) and (4): <epsilon> = (Vm*I + Vf*A)*<epsilon_m>
    #    -> <epsilon_m> = inv(Vm*I + Vf*A)*<epsilon>
    # 6. From (1), (3), and (4): <sigma> = (Vm*Cm + Vf*Cf*A)*<epsilon_m>
    # 7. Substitute <epsilon_m> from (5) into (6):
    #    <sigma> = (Vm*Cm + Vf*Cf*A)*inv(Vm*I + Vf*A)*<epsilon>
    # 8. Since <sigma> = C*<epsilon>, we get the final expression for C.

    print("The expression for the effective average elastic moduli (C) is:")
    # We use '*' to denote tensor contraction (e.g., A:B or A_ijkl*B_klmn)
    # and inv() to denote the inverse of a fourth-order tensor.
    final_equation = f"{C} = ({Vm} * {Cm} + {Vf} * {Cf} * {A}) * inv({Vm} * {I} + {Vf} * {A})"
    print(final_equation)
    print("\nWhere the terms in the equation are:")
    print(f"{C}: The effective average elastic moduli of the composite")
    print(f"{Cm}: The fourth-order elasticity tensor of the matrix")
    print(f"{Cf}: The fourth-order elasticity tensor of the fiber")
    print(f"{Vm}: The volume fraction of the matrix")
    print(f"{Vf}: The volume fraction of the fiber")
    print(f"{A}: The Eshelby strain-concentration tensor")
    print(f"{I}: The fourth-order identity tensor")
    print("inv(...): Represents the inverse of the fourth-order tensor in the parenthesis")

mori_tanaka_equation()
<<<C = (Vm * Cm + Vf * Cf * A) * inv(Vm * I + Vf * A)>>>