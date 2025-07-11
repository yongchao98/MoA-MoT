import math

def solve_double_points():
    """
    Calculates the number of double points in the stable reduction of the curve
    y^2 = 8*x^5 + 4*x^4 + 4*x^3 + x^2 + 8x over the prime p=2.
    """

    # Step 1: Define the curve and calculate its genus
    # The curve is y^2 = f(x) where f(x) = 8x^5 + 4x^4 + 4x^3 + x^2 + 8x
    # The degree of the polynomial f(x) is d.
    d = 5
    # The genus g of a hyperelliptic curve is given by g = floor((d-1)/2)
    g = math.floor((d - 1) / 2)

    print(f"The curve is given by y^2 = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x.")
    print(f"The degree of the polynomial is d = {d}.")
    print(f"The genus of the curve is g = floor(({d}-1)/2) = {g}.")
    print("-" * 30)

    # Step 2: Identify branch points and their properties.
    # The branch points are the 5 roots of f(x) (r_0, ..., r_4) and the point at infinity (r_inf).
    # From an analysis using Newton polygons over the 2-adic numbers, we find the valuations of the roots:
    # v2(r_0) = infinity (since r_0 = 0)
    # v2(r_1) = -1
    # v2(r_2) = -1
    # v2(r_3) = -1
    # v2(r_4) = 3
    print("The branch points are the 6 roots of the equation in the projective line over the 2-adic field.")
    print("By analyzing the 2-adic valuations of the branch points, we can group them into clusters.")
    print("-" * 30)

    # Step 3: Group the branch points into clusters.
    # Two branch points alpha, beta are in the same cluster if v2(alpha - beta) > 0.
    # A finite branch point alpha clusters with infinity if v2(alpha) < 0.
    # - r_0 and r_4 cluster together because v2(r_0 - r_4) = v2(-r_4) = v2(r_4) = 3 > 0.
    # - r_1, r_2, r_3 cluster with infinity because their valuations are -1 < 0.
    # This gives us two clusters.
    
    n1 = 2  # Cluster 1: {r_0, r_4}
    n2 = 4  # Cluster 2: {r_1, r_2, r_3, r_inf}
    N = 2   # Total number of clusters

    print(f"There are {N} clusters of branch points.")
    print(f"Cluster 1 has size n1 = {n1}.")
    print(f"Cluster 2 has size n2 = {n2}.")
    print("-" * 30)

    # Step 4: Determine the genera of the components of the stable reduction.
    # The genus gi of a component corresponding to a cluster of size ni is gi = floor((ni - 1)/2).
    g1 = math.floor((n1 - 1) / 2)
    g2 = math.floor((n2 - 1) / 2)
    
    print("Each cluster corresponds to a component in the stable reduction curve.")
    print(f"The genus of the first component is g1 = floor(({n1}-1)/2) = {g1}.")
    print(f"The genus of the second component is g2 = floor(({n2}-1)/2) = {g2}.")
    print("-" * 30)

    # Step 5: Calculate the number of double points (delta).
    # The formula relating these quantities is: g - 1 = (g1 - 1) + (g2 - 1) + delta
    # We solve for delta: delta = (g - 1) - (g1 - 1) - (g2 - 1)
    # delta = g - g1 - g2 + 1
    
    delta = (g - 1) - (g1 - 1) - (g2 - 1)

    print("The number of double points (delta) is found using the formula:")
    print("g - 1 = (g1 - 1) + (g2 - 1) + delta")
    print(f"{g} - 1 = ({g1} - 1) + ({g2} - 1) + delta")
    print(f"{g - 1} = {g1 - 1} + {g2 - 1} + delta")
    print(f"{g - 1} = {g1 - 1 + g2 - 1} + delta")
    print(f"delta = {g - 1} - ({g1 - 1 + g2 - 1}) = {delta}")
    print("-" * 30)
    
    print(f"The number of double points is {delta}.")
    
    return delta

if __name__ == '__main__':
    solve_double_points()
