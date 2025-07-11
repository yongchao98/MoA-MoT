def solve_mori_tanaka_expression():
    """
    This function constructs and prints the expression for the effective
    elastic moduli (C) using the Mori-Tanaka model.
    """
    # Define the symbols used in the model
    C = "C"
    Cm = "Cm"
    Cf = "Cf"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"
    I = "I"

    # The relationship between volume fractions is Vm + Vf = 1
    print(f"The relationship between the volume fractions is: {Vm} + {Vf} = 1")
    print("-" * 50)

    # The final expression for the effective average elastic moduli C
    # It is commonly written as: C = Cm + Vf * (Cf - Cm) * A * (Vm*I + Vf*A)^-1
    # where the term A * (Vm*I + Vf*A)^-1 is the strain concentration tensor
    # that relates the average fiber strain to the average composite strain.

    # We will print each "number" or symbol in the final equation as requested.
    print("The Mori-Tanaka expression for the effective average elastic moduli (C) is composed of the following terms:")

    part1 = Cm
    part2 = Vf
    part3 = f"({Cf} - {Cm})"
    part4 = A
    part5 = f"({Vm}*{I} + {Vf}*{A})**-1" # Using **-1 to denote the tensor inverse

    print(f"Final Equation: {C} = {part1} + {part2} * {part3} * {part4} * {part5}")
    print("\nWhere:")
    print(f"* {C}: The effective average elastic moduli tensor of the composite.")
    print(f"* {Cm}: The elasticity tensor of the matrix.")
    print(f"* {Vf}: The volume fraction of the fibers.")
    print(f"* {Cf}: The elasticity tensor of the fiber.")
    print(f"* {A}: The Eshelby strain-concentration tensor.")
    print(f"* {Vm}: The volume fraction of the matrix.")
    print(f"* {I}: The fourth-order identity tensor.")
    print(f"* (**-1): Denotes the inverse of the tensor.")


solve_mori_tanaka_expression()

final_answer = "C = Cm + Vf * (Cf - Cm) * A * (Vm*I + Vf*A)**-1"
print(f"\n<<<C = Cm + Vf * (Cf - Cm) * A * (Vm*I + Vf*A)**-1>>>")