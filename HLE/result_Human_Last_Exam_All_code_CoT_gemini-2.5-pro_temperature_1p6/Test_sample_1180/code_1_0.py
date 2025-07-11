def solve():
    """
    This script explains the steps to calculate the thickness of the double point.
    """
    
    # Genus of the original curve z^2 = 2x^5 + 2x^3 + 1.
    # The degree of the polynomial in x is 5.
    # Genus g = floor((deg-1)/2)
    g = 2

    # After a change of variables (x=t/2, w=4z), the curve is birationally
    # equivalent to w^2 = t^5 + 4t^3 + 16.
    # The Newton Polygon of t^5 + 4t^3 + 16 shows the roots of the polynomial
    # split into two clusters.
    # Cluster 1 (R1): 3 roots.
    # Cluster 2 (R2): 2 roots.
    num_roots_R1 = 3
    num_roots_R2 = 2

    # The stable reduction has two main components, C1 and C2, corresponding
    # to these clusters of roots.
    # The genus of the component for a cluster of k roots is g_k = floor((k-1)/2).
    g1 = (num_roots_R1 - 1) // 2  # Genus of component C1
    g2 = (num_roots_R2 - 1) // 2  # Genus of component C2

    # The arithmetic genus of the stable reduction fiber (pa) must equal the genus
    # of the original curve (g).
    pa = g

    # The formula for the arithmetic genus of a curve degenerating into two components is:
    # pa = g1 + g2 + thickness - 1
    # We solve for the thickness.
    # thickness = pa - g1 - g2 + 1
    thickness = pa - g1 - g2 + 1
    
    print("The genus of the curve is g = {}.".format(g))
    print("The stable reduction has two components, C1 and C2.")
    print("The number of roots for C1 is {}.".format(num_roots_R1))
    print("The number of roots for C2 is {}.".format(num_roots_R2))
    print("The genus of C1 is g1 = {}.".format(g1))
    print("The genus of C2 is g2 = {}.".format(g2))
    print("The arithmetic genus of the stable fiber is pa = {}.".format(pa))
    print("The relationship is given by the formula: pa = g1 + g2 + thickness - 1")
    print("Substituting the values:")
    print("{} = {} + {} + thickness - {}".format(pa, g1, g2, 1))
    print("Solving for the thickness gives: thickness = {}".format(thickness))

solve()