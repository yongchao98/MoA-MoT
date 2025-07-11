def solve():
    """
    Calculates the number of double points in the stable reduction of the curve
    y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5 above 2.
    """

    # Step 1: Genus of the generic curve y^2 = f(x)
    # The degree of f(x) is 5.
    d = 5
    # The genus of the generic curve is g = floor((d-1)/2).
    g = (d - 1) // 2
    
    # Step 2 & 3: Find a model with a smooth reduction and calculate its genus.
    # The process of finding a stable reduction involves a sequence of transformations.
    # As described in the explanation, after 3 transformations, we arrive at a model
    # whose reduction modulo 2 is y^2 = x^2 + x.
    # This reduced curve is smooth and rational (isomorphic to a parabola).
    # The genus of a smooth rational curve is 0.
    g_s = 0

    # Step 4: The number of double points (delta) in the stable reduction is
    # given by the genus drop formula: delta = g - g_s.
    delta = g - g_s

    # The final equation is delta = g - g_s
    print(f"The genus of the original curve is g = {g}.")
    print(f"The genus of the stable reduction's component is g_s = {g_s}.")
    print(f"The number of double points is delta = g - g_s = {g} - {g_s} = {delta}.")

solve()