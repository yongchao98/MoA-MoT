def calculate_symmetry_breaking():
    """
    Calculates the number of massive gauge bosons resulting from the
    spontaneous symmetry breaking of SU(3) -> SU(2) x U(1).
    """

    # 1. Calculate generators for the initial group G = SU(3)
    n_g = 3
    dim_g = n_g**2 - 1
    print(f"Initial group G = SU({n_g})")
    print(f"Number of generators for G: {n_g}^2 - 1 = {dim_g}\n")

    # 2. Calculate generators for the residual group H = SU(2) x U(1)
    n_h_su2 = 2
    dim_h_su2 = n_h_su2**2 - 1
    dim_h_u1 = 1  # U(1) has one generator
    dim_h = dim_h_su2 + dim_h_u1
    print("Residual group H = SU(2) x U(1)")
    print(f"Number of generators for SU(2): {n_h_su2}^2 - 1 = {dim_h_su2}")
    print(f"Number of generators for U(1): {dim_h_u1}")
    print(f"Total number of unbroken generators (dim(H)): {dim_h_su2} + {dim_h_u1} = {dim_h}\n")

    # 3. Calculate the number of broken generators
    broken_generators = dim_g - dim_h
    print("The number of broken generators is dim(G) - dim(H).")
    print(f"Number of broken generators = {dim_g} - {dim_h} = {broken_generators}\n")

    # 4. Relate broken generators to massive gauge bosons
    print("In a gauge theory, each broken generator corresponds to a gauge boson")
    print("that acquires mass via the Higgs mechanism.")
    print(f"Therefore, there are {broken_generators} massive gauge bosons.")

    # Final equation breakdown
    print("\n--- Final Equation ---")
    print("Number of Massive Gauge Bosons = (Generators of SU(3)) - (Generators of SU(2) + Generators of U(1))")
    print(f"                               = ({dim_g}) - ({dim_h_su2} + {dim_h_u1})")
    print(f"                               = {dim_g} - {dim_h}")
    print(f"                               = {broken_generators}")

calculate_symmetry_breaking()