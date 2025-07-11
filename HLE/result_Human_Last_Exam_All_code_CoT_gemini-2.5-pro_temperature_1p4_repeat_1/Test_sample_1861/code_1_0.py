def solve_manifold_problem():
    """
    Solves the differential geometry problem by presenting a step-by-step logical argument.
    """
    print("Problem Analysis:")
    print("Let M be a 2-dim orientable manifold (torus, cylinder, or plane).")
    print("Let η be a 1-form on M.")
    print("The key condition is: for any x, y in M, there is a diffeomorphism F with F(x) = y and F*η = η.")
    print("This means the group G of η-preserving diffeomorphisms acts transitively on M.")
    print("-" * 30)

    print("Step 1: Determine the nature of dη.")
    print("Since F*η = η for any F in G, we can apply the exterior derivative d:")
    print("d(F*η) = d(η)")
    print("Because d commutes with F*, this gives F*(dη) = dη.")
    print("So, the 2-form dη is also invariant under the group G.")
    print("Since G acts transitively, dη must be a constant 2-form.")
    print("We can write dη = C * ω, where C is a constant and ω is a G-invariant volume form.")
    print("The main question is: Must C be zero?")
    print("-" * 30)

    def analyze_torus():
        """Analysis for the 2-torus."""
        print("Step 2: Analysis for M = 2-torus (T²)")
        print("The torus is a compact manifold with no boundary (∂T² = ∅).")
        print("By Stokes' Theorem: ∫_M dα = ∫_∂M α.")
        print("Applying this to η on the torus M = T²:")
        print("∫_T² dη = ∫_(∂T²) η = 0.")
        print("We also have ∫_T² dη = ∫_T² (C * ω) = C * Area(T²).")
        print("So, C * Area(T²) = 0.")
        print("Since Area(T²) is a finite positive number, the constant C must be 0.")
        print("Conclusion for Torus: It is necessary that C = 0, which means dη = 0.")
        print("-" * 30)

    def analyze_cylinder():
        """Analysis for the cylinder."""
        print("Step 3: Analysis for M = Cylinder (S¹ x ℝ)")
        print("The cylinder is non-compact. A transitive group G must contain translations along the ℝ-axis.")
        print("This means there's a vector field v (like ∂/∂t) whose flow preserves η, so the Lie derivative L_v(η) = 0.")
        print("Using Cartan's formula, L_v(η) = i_v(dη) + d(i_v(η)), we get that i_v(dη) must be an exact 1-form.")
        print("Let dη = C * dθ ∧ dt. The vector field for ℝ-translation is v = ∂/∂t.")
        print("i_v(dη) = i_{∂/∂t}(C * dθ ∧ dt) = C * dθ.")
        print("The 1-form dθ is closed on the cylinder but NOT exact. It generates the cohomology group H¹(S¹ x ℝ).")
        print("For the form 'C * dθ' to be exact, the constant C must be 0.")
        print("Conclusion for Cylinder: It is necessary that C = 0, which means dη = 0.")
        print("-" * 30)

    def analyze_plane():
        """Analysis for the plane."""
        print("Step 4: Analysis for M = Plane (ℝ²)")
        print("The plane is non-compact. A transitive group G on ℝ² must contain translations.")
        print("This means η must be invariant under vector fields v_x = ∂/∂x and v_y = ∂/∂y.")
        print("The condition L_{v_x}(η) = 0 and L_{v_y}(η) = 0 means η is invariant under all translations.")
        print("Let η = A(x,y)dx + B(x,y)dy.")
        print("Invariance under x-translations implies A and B are functions of y only.")
        print("Invariance under y-translations implies A and B are functions of x only.")
        print("For both to be true, A and B must be constants.")
        print("So, η = A*dx + B*dy, where A and B are constants.")
        print("The exterior derivative is dη = d(A*dx + B*dy) = 0.")
        print("Conclusion for Plane: It is necessary that dη = 0.")
        print("-" * 30)

    analyze_torus()
    analyze_cylinder()
    analyze_plane()

    print("Final Conclusion:")
    print("In all three cases (torus, cylinder, and plane), the given condition necessarily implies that dη must be zero.")
    final_number = 0
    print(f"The final equation is: dη = {final_number}")

solve_manifold_problem()