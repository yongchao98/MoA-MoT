def calculate_genus_of_pentagon_linkage():
    """
    This function calculates the genus of the configuration space of a hinged
    regular pentagon with two adjacent vertices nailed to the plane.

    The configuration space is a smooth surface whose topology is known.
    Its genus 'g' can be determined from its first Betti number, b1, which
    represents the number of independent, non-trivial loops on the surface.
    The relationship for a connected, orientable surface is g = b1 / 2.

    From established results in algebraic topology (Kamiyama, 2000), the Betti
    numbers for this specific configuration space are:
    b0 = 1 (it's one connected component)
    b1 = 8
    b2 = 1 (it's a closed, orientable surface)
    """

    # The first Betti number of the configuration space.
    b1 = 8

    # The formula relating the genus (g) to the first Betti number (b1).
    # g = b1 / 2
    genus = b1 / 2

    print("The problem is to find the genus of the configuration space of a specific pentagonal linkage.")
    print("This space is a smooth, orientable surface.")
    print("The genus 'g' is related to the first Betti number 'b1' by the formula: g = b1 / 2.")
    print(f"The first Betti number for this surface is known to be: b1 = {b1}")
    print("Calculating the genus:")
    print(f"g = {b1} / 2 = {int(genus)}")
    print(f"The genus of the surface is {int(genus)}.")

calculate_genus_of_pentagon_linkage()