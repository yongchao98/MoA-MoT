def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space X.

    The space X is the disjoint union of:
    1. The torus (T^2)
    2. The sphere (S^2)
    3. The real line (R)
    4. A three-point discrete space (D_3)
    5. A five-point discrete space (D_5)

    The equivalence relation x ~ y holds if there is a homeomorphism f: X -> X
    such that f(x) = y.
    """
    
    # A homeomorphism on X must map its components to homeomorphic components.
    # We first check if any of the 5 components are homeomorphic to each other.
    # - Torus vs Sphere: No, different fundamental groups.
    # - Compact spaces (Torus, Sphere) vs Non-compact (Real line): No.
    # - Connected spaces (Torus, Sphere, Real line) vs Discrete (not connected): No.
    # - 3-point discrete vs 5-point discrete: No, different number of points.
    # Conclusion: No two components are homeomorphic. Thus, any homeomorphism of X
    # must map each component to itself. This means equivalence classes are contained
    # within single components.

    # We now count the number of equivalence classes within each component.
    # This is equivalent to checking if each component space is homogeneous.
    
    # 1. The Torus: It is a homogeneous space. Any point can be mapped to any
    #    other via a translation, which is a homeomorphism.
    classes_torus = 1
    
    # 2. The Sphere: It is a homogeneous space. Any point can be mapped to any
    #    other via a rotation, which is a homeomorphism.
    classes_sphere = 1
    
    # 3. The Real Line: It is a homogeneous space. Any point can be mapped to
    #    any other via a translation (e.g., f(z) = z + (y-x)), which is a homeomorphism.
    classes_real_line = 1
    
    # 4. A three-point discrete space: Any permutation of points is a homeomorphism.
    #    The group of permutations acts transitively, so the space is homogeneous.
    classes_d3 = 1
    
    # 5. A five-point discrete space: Similar to the 3-point space, it's homogeneous.
    classes_d5 = 1

    # The total number of equivalence classes in X is the sum of the number of
    # classes in each of its components.
    total_classes = (classes_torus + 
                     classes_sphere + 
                     classes_real_line + 
                     classes_d3 + 
                     classes_d5)
                     
    print("Number of equivalence classes for each component:")
    print(f"- Torus: {classes_torus}")
    print(f"- Sphere: {classes_sphere}")
    print(f"- Real line: {classes_real_line}")
    print(f"- Three-point discrete space: {classes_d3}")
    print(f"- Five-point discrete space: {classes_d5}")
    
    print("\nThe total number of equivalence classes is the sum:")
    print(f"{classes_torus} + {classes_sphere} + {classes_real_line} + {classes_d3} + {classes_d5} = {total_classes}")

solve_equivalence_classes()