def calculate_symmetry_breaking():
    """
    Calculates the number of broken generators for SU(3) -> SU(2) x U(1)
    and determines the number of resulting massive gauge bosons.
    """
    # Number of generators for the initial group G = SU(3)
    n_g = 3
    generators_g = n_g**2 - 1

    # Number of generators for the residual group H = SU(2) x U(1)
    # For SU(2)
    n_h_su2 = 2
    generators_h_su2 = n_h_su2**2 - 1
    # For U(1)
    generators_h_u1 = 1
    # Total for H
    generators_h = generators_h_su2 + generators_h_u1

    # Number of broken generators
    broken_generators = generators_g - generators_h

    print("In the context of a non-Abelian gauge theory, spontaneous symmetry breaking")
    print("leads to massive gauge bosons corresponding to the broken generators.")
    print("\nStep 1: Calculate generators of the initial group G = SU(3)")
    print(f"dim(SU(3)) = {n_g}^2 - 1 = {generators_g}")

    print("\nStep 2: Calculate generators of the residual group H = SU(2) x U(1)")
    print(f"dim(SU(2)) = {n_h_su2}^2 - 1 = {generators_h_su2}")
    print(f"dim(U(1)) = {generators_h_u1}")
    print(f"dim(H) = {generators_h_su2} + {generators_h_u1} = {generators_h}")

    print("\nStep 3: Calculate the number of broken generators (dim(G) - dim(H))")
    print("This number corresponds to the number of massive gauge bosons.")
    print(f"Final Equation: {generators_g} - ({generators_h_su2} + {generators_h_u1}) = {broken_generators}")

    print(f"\nConclusion: The breaking SU(3) -> SU(2) x U(1) results in {broken_generators} massive gauge bosons.")

calculate_symmetry_breaking()