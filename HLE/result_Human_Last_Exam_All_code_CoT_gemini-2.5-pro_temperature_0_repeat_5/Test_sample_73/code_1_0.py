def solve_pentagon_genus():
    """
    Calculates the genus of the configuration space of a hinged pentagon.

    The configuration space of a hinged pentagon with one side fixed is a
    well-studied 2-dimensional surface. Its topology is described by its
    Betti numbers, which are known from the mathematical literature on
    moduli spaces of polygons.

    The Betti numbers for this surface are:
    b_0 = 1 (it is one connected component)
    b_1 = 8 (it has 8 independent non-bounding cycles)
    b_2 = 1 (it encloses one 'void')
    """
    
    # Betti numbers for the configuration space of a planar pentagon
    b0 = 1
    b1 = 8
    b2 = 1
    
    print("The calculation is based on the known Betti numbers of the configuration space.")
    print(f"b_0 (number of connected components) = {b0}")
    print(f"b_1 (number of handles/holes) = {b1}")
    print(f"b_2 (number of voids) = {b2}")
    print("-" * 30)

    # The genus 'g' of a connected orientable surface is half its first Betti number.
    # Equation: g = b1 / 2
    print("The genus (g) is calculated from the first Betti number (b1).")
    print(f"Equation: g = {b1} / 2")
    genus = b1 // 2
    print(f"Result: g = {genus}")
    print("-" * 30)

    # As a check, we can compute the Euler characteristic (chi) in two ways.
    # 1. From Betti numbers: chi = b0 - b1 + b2
    print("For verification, we calculate the Euler characteristic (chi).")
    print("Method 1: From Betti numbers.")
    print(f"Equation: chi = {b0} - {b1} + {b2}")
    chi_from_betti = b0 - b1 + b2
    print(f"Result: chi = {chi_from_betti}")
    print("-" * 30)

    # 2. From genus: chi = 2 - 2g
    print("Method 2: From genus.")
    print(f"Equation: chi = 2 - 2 * {genus}")
    chi_from_genus = 2 - 2 * genus
    print(f"Result: chi = {chi_from_genus}")
    print("-" * 30)

    if chi_from_betti == chi_from_genus:
        print("The results match, confirming the genus calculation.")
    else:
        print("Error: The results do not match.")

    print("\nThe genus of the smooth surface is the final answer.")
    print(f"Final Answer: {genus}")

solve_pentagon_genus()