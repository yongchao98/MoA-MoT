def solve_manifold_identity():
    """
    This function simplifies the group presentation derived from the Heegaard diagram
    to identify the fundamental group and thus the manifold.
    """
    print("Step 1: The fundamental group presentation")
    print("The diagram is interpreted to give the following presentation for the fundamental group π₁(M):")
    print("Generators: x₁, x₂, x₃")
    print("Relators:")
    print("  r₁: x₂ * x₃ = 1")
    print("  r₂: x₁ * x₃ = 1")
    print("  r₃: x₁ * x₂ = 1")
    print("-" * 40)

    print("Step 2: Simplify the relations")
    print("From r₁, x₂x₃ = 1, we can solve for x₃:")
    print("  x₃ = x₂⁻¹")
    print("\nFrom r₂, x₁x₃ = 1, we can also solve for x₃:")
    print("  x₃ = x₁⁻¹")
    print("\nEquating the two expressions for x₃ gives:")
    print("  x₂⁻¹ = x₁⁻¹")
    print("Taking the inverse of both sides, we find:")
    print("  x₂ = x₁")
    print("-" * 40)

    print("Step 3: Substitute and find the final relation")
    print("Substitute x₂ = x₁ into the third relation, r₃: x₁x₂ = 1")
    print("  x₁ * (x₁) = 1")
    print("This simplifies to the final relation:")
    final_equation_lhs = "x₁²"
    final_equation_rhs = "1"
    print(f"  {final_equation_lhs} = {final_equation_rhs}")
    print("-" * 40)

    print("Step 4: Identify the group and the manifold")
    print("The entire group can be described by a single generator x₁ with the relation x₁² = 1.")
    print("The other generators are x₂ = x₁ and x₃ = x₁⁻¹ = x₁.")
    print("The fundamental group is π₁(M) = <x₁ | x₁² = 1>, which is the cyclic group of order 2, Z₂.")
    print("\nA 3-manifold with this fundamental group is the lens space L(2,1).")
    print("This manifold is commonly known as real projective 3-space, or RP³.")

solve_manifold_identity()