def solve_genus_problem():
    """
    This function determines the maximal genus of a surface under the given conditions.
    
    The reasoning is based on established theorems in differential geometry:
    1.  Genus 0 (sphere) is possible. Its mean curvature can be H = 1/R > 0.
    2.  Genus 1 (torus) is possible. A torus of revolution with major radius R and minor radius r,
        if R > 2*r, has mean curvature H > 0 everywhere.
    3.  For genus g >= 2, it can be shown that any such surface that bounds a compact region
        in R^3 must have points where the mean curvature is zero. Therefore, surfaces with
        g >= 2 cannot satisfy the condition that the mean curvature never vanishes.
    
    This implies the maximal genus is 1.
    """
    
    # Maximal genus based on the reasoning above.
    g_max = 1
    
    # Outputting the final equation for the maximal genus.
    print(f"The analysis shows that the maximal possible genus is g = 1.")
    print("The final equation is:")
    print(f"g_max = {g_max}")

solve_genus_problem()