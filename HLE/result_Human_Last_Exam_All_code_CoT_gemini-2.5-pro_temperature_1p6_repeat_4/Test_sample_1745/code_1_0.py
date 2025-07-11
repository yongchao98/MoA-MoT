def solve_edl_potential():
    """
    This function derives and prints the expression for the Electrical Double-Layer (EDL)
    potential distribution in a parallel-plate microchannel under the specified conditions.
    """

    # Define string representations of the variables for clear printing.
    psi = "ψ(y)"
    z_a1 = "z_a1"
    z_1 = "z₁"
    beta = "β"
    k = "k"
    H = "H"
    y = "y"

    # Print the step-by-step derivation
    print("Derivation of the Electrical Double-Layer (EDL) Potential Distribution")
    print("=" * 75)

    print("\nStep 1: The Governing Equation")
    print("The EDL potential distribution, ψ(y), is governed by the linearized Poisson-Boltzmann equation:")
    print(f"  d²ψ/dy² = {k}² * ψ")
    print("where 'k' is the Debye-Hückel parameter.")

    print("\nStep 2: Boundary Conditions")
    print("The microchannel is defined by two parallel plates at y=0 and y=H.")
    print(f"- At the bottom plate (y=0): ψ(0) = {z_a1} = {z_1} * (1 + {beta}*{k})")
    print(f"- At the top plate (y=H): The zeta potential z₂ is 0, so ψ(H) = 0")

    print("\nStep 3: General Solution")
    print("A convenient form of the general solution to this ODE is:")
    print(f"  {psi} = A * sinh({k}*({H}-{y})) + B * cosh({k}*({H}-{y}))")
    print("where A and B are constants to be determined.")

    print("\nStep 4: Applying Boundary Conditions")
    print("First, we apply the condition at the top plate, ψ(H) = 0:")
    print(f"  0 = A * sinh({k}*({H}-{H})) + B * cosh({k}*({H}-{H}))")
    print(f"  0 = A * sinh(0) + B * cosh(0)")
    print(f"  0 = A * 0 + B * 1   =>   This gives B = 0")
    print(f"\nThe solution simplifies to: {psi} = A * sinh({k}*({H}-{y}))")

    print("\nNext, we apply the condition at the bottom plate, ψ(0) = z_a1:")
    print(f"  {z_a1} = A * sinh({k}*({H}-0))")
    print(f"  {z_1}*(1 + {beta}*{k}) = A * sinh({k}*{H})")
    print(f"  Solving for A: A = ({z_1}*(1 + {beta}*{k})) / sinh({k}*{H})")

    print("\nStep 5: Final Expression")
    print("Substituting the expression for A back into the simplified solution yields the final result.")
    print("The final equation for the EDL potential distribution is:")
    print("-" * 60)
    print(f"        {z_1} * (1 + {beta}*{k}) * sinh({k}*({H} - {y}))")
    print(f"{psi} = -------------------------------------------------")
    print(f"                      sinh({k}*{H})")
    print("-" * 60)

    print("\nComponents of the Final Equation:")
    print(f"  {psi}: The electrical potential at a transverse position y.")
    print(f"  {z_1}:   The base zeta potential at the bottom surface (y=0).")
    print(f"  {beta}:    The slip length.")
    print(f"  {k}:     The Debye–Huckel parameter.")
    print(f"  {H}:     The height of the microchannel.")
    print(f"  {y}:     The transverse coordinate, measured from the bottom plate.")
    print("  sinh:  The hyperbolic sine function.")

if __name__ == '__main__':
    solve_edl_potential()