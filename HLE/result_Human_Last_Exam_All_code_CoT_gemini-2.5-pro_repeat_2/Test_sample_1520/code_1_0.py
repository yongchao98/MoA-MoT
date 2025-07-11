def solve_symmetry_breaking():
    """
    Calculates the number of massive gauge bosons resulting from the spontaneous
    symmetry breaking of SU(3) to SU(2) x U(1).
    """

    # Step 1: Calculate generators for the initial group G = SU(3)
    N_initial = 3
    generators_G = N_initial**2 - 1
    print(f"The initial symmetry group is G = SU({N_initial}).")
    print(f"Number of generators for SU({N_initial}) = {N_initial}^2 - 1 = {generators_G}.")
    print("-" * 30)

    # Step 2: Calculate generators for the residual group H = SU(2) x U(1)
    N_residual_su2 = 2
    generators_su2 = N_residual_su2**2 - 1
    generators_u1 = 1
    generators_H = generators_su2 + generators_u1
    print(f"The residual symmetry group is H = SU({N_residual_su2}) x U(1).")
    print(f"Number of generators for SU({N_residual_su2}) = {N_residual_su2}^2 - 1 = {generators_su2}.")
    print(f"Number of generators for U(1) = {generators_u1}.")
    print(f"Total number of unbroken generators for H = {generators_su2} + {generators_u1} = {generators_H}.")
    print("-" * 30)

    # Step 3 & 4: Calculate the number of broken generators, which corresponds
    # to the number of massive gauge bosons in a gauge theory.
    broken_generators = generators_G - generators_H
    print("In a gauge theory, the number of massive gauge bosons is equal to the number of broken generators.")
    print("Number of broken generators = (Generators of G) - (Generators of H)")
    print(f"The final calculation is: {generators_G} - ({generators_su2} + {generators_u1}) = {broken_generators}")
    print("-" * 30)
    print(f"Result: There are {broken_generators} massive gauge bosons.")
    print("This corresponds to option E.")

solve_symmetry_breaking()