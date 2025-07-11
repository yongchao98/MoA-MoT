def solve_heegaard_diagram():
    """
    Analyzes the given Heegaard diagram to identify the corresponding 3-manifold.
    """
    print("This program determines the 3-manifold represented by the Heegaard diagram.")
    print("-" * 50)

    # Step 1: Analyze the Heegaard diagram's structure.
    genus = 3
    print(f"Step 1: Analyzing the Heegaard Diagram")
    print(f"The diagram is a genus-{genus} Heegaard diagram. The red curves α₁, α₂, α₃ represent")
    print("a basis for the homology of one handlebody, while the blue curves β₁, β₂, β₃")
    print("represent the attaching curves for the second handlebody.")
    print("-" * 50)

    # Step 2: Set up the fundamental group presentation.
    print("Step 2: Deriving the Fundamental Group π₁(M)")
    print("The fundamental group is given by a presentation with one generator for each α-curve")
    print("and one relator for each β-curve.")
    print(f"Generators: x₁, x₂, x₃")
    print(f"Relators: r₁, r₂, r₃")
    print("The presentation is: π₁(M) = <x₁, x₂, x₃ | r₁, r₂, r₃>")
    print("-" * 50)
    
    # Step 3: Determine the relators by tracing the β-curves.
    print("Step 3: Reading the Relators from the Diagram")
    print("We trace each β-curve and record its intersections with the α-curves.")
    # Relator 1 from β₁
    relator1_eq = "x₃ * x₂⁻¹ = 1"
    relator1_lhs = "x₃"
    relator1_rhs = "x₂"
    print(f" - The curve β₁ gives the relator r₁: {relator1_eq}, which simplifies to {relator1_lhs} = {relator1_rhs}.")

    # Relator 2 from β₂
    relator2_eq = "x₁ * x₃⁻¹ = 1"
    relator2_lhs = "x₁"
    relator2_rhs = "x₃"
    print(f" - The curve β₂ gives the relator r₂: {relator2_eq}, which simplifies to {relator2_lhs} = {relator2_rhs}.")

    # Relator 3 from β₃
    relator3_eq = "x₂ * x₁⁻¹ = 1"
    relator3_lhs = "x₂"
    relator3_rhs = "x₁"
    print(f" - The curve β₃ gives the relator r₃: {relator3_eq}, which simplifies to {relator3_lhs} = {relator3_rhs}.")
    print("-" * 50)

    # Step 4: Simplify the group presentation.
    print("Step 4: Simplifying the Group Presentation")
    final_relation_var1 = "x₁"
    final_relation_var2 = "x₂"
    final_relation_var3 = "x₃"
    print("Combining the relations, we get:")
    print(f"  {final_relation_var1} = {final_relation_var2} = {final_relation_var3}")
    print("This means the group has a single generator 'x' with no relations.")
    print("The fundamental group is π₁(M) = <x> ≅ ℤ (the infinite cyclic group).")
    print("-" * 50)

    # Step 5: Identify the manifold based on its fundamental group.
    print("Step 5: Identifying the 3-Manifold")
    print("A classification theorem in topology states that a closed, orientable 3-manifold")
    print("with fundamental group ℤ is homeomorphic to S² × S¹.")
    print("-" * 50)
    
    # Final Conclusion
    manifold_name = "S² × S¹"
    manifold_description = "the product of a 2-sphere and a circle"
    print(f"Conclusion: The Heegaard diagram represents the three-manifold {manifold_name},")
    print(f"which is {manifold_description}.")

if __name__ == "__main__":
    solve_heegaard_diagram()