def solve_phi4_fixed_point():
    """
    Calculates and explains the leading order expression for the Wilson-Fisher
    fixed point coupling in phi^4 theory near d=4 dimensions.
    """
    print("Derivation of the fixed point coupling u* in ϕ⁴ theory.")
    print("============================================================")

    print("\n1. The Context")
    print("   We analyze the ϕ⁴ theory using the Renormalization Group (RG) in")
    print("   d = 4 - ϵ dimensions, where ϵ is a small parameter.")
    print("   The behavior of the theory under changes in energy scale is governed")
    print("   by the beta function, β(u), for the dimensionless coupling u.")

    print("\n2. The Beta Function")
    print("   To one-loop order, the beta function for the coupling u is:")
    print("   β(u) = -ϵ*u + (3 / (16*π²)) * u²")
    print("   - The '-ϵ*u' term arises from the canonical scaling of the coupling.")
    print("   - The '(3 / (16*π²)) * u²' term is the quantum correction from one-loop diagrams.")

    print("\n3. Finding the Fixed Point")
    print("   A fixed point, denoted u*, is a point where the coupling does not change")
    print("   with the energy scale, meaning the beta function is zero.")
    print("   We solve the equation: β(u*) = 0")
    print("   -ϵ*u* + (3 / (16*π²)) * (u*)² = 0")

    print("\n4. Solving the Equation")
    print("   We can factor the equation as:")
    print("   u* * [-ϵ + (3 / (16*π²)) * u*] = 0")
    print("   This gives two solutions:")
    print("   a) u* = 0 (The Gaussian or trivial fixed point)")
    print("   b) -ϵ + (3 / (16*π²)) * u* = 0 (The Wilson-Fisher or non-trivial fixed point)")
    print("\n   Solving for the non-trivial fixed point:")
    print("   (3 / (16*π²)) * u* = ϵ")
    print("   u* = ϵ * (16*π² / 3)")

    print("\n5. Final Expression")
    print("   The leading order expression for the fixed point coupling u* is:")

    # As requested, output each number in the final equation
    num_16 = 16
    exponent = 2
    num_3 = 3
    print(f"\n   u* = ({num_16} * π^{exponent} / {num_3}) * ϵ\n")
    print("============================================================")

if __name__ == '__main__':
    solve_phi4_fixed_point()