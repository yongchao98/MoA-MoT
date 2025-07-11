import math

def solve_curve_reduction():
    """
    Calculates the number of double points in the stable reduction of the curve
    y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5 above p=2.
    """
    
    # Step 1: Define the curve and find its genus
    # The polynomial is f(x) = 8x^5 + 4x^4 + 4x^3 + x^2 + 8x
    # The degree of f(x) is 5.
    deg_f = 5
    g = math.floor((deg_f - 1) / 2)
    print(f"The curve is y^2 = 8x^5 + 4x^4 + 4x^3 + x^2 + 8x.")
    print(f"The genus of the curve is g = floor(({deg_f} - 1) / 2) = {g}.")
    print("-" * 20)

    # Step 2 & 3: Explain the method based on branch points
    print("To find the number of double points, we analyze the 2-adic valuations of the curve's branch points.")
    print("The branch points are the roots of f(x) and the point at infinity (since the degree is odd).")
    print("-" * 20)

    # Step 4: Find the valuations of the roots of f(x)
    # f(x) = x * (8x^4 + 4x^3 + 4x^2 + x + 8)
    # One root is alpha_1 = 0. Its 2-adic valuation is +infinity.
    print("One root is x = 0.")
    # For the other roots, we use the Newton Polygon of g(x) = 8x^4 + 4x^3 + 4x^2 + x + 8 at p=2.
    # The coefficients of g(x) are d_4=8, d_3=4, d_2=4, d_1=1, d_0=8.
    # The 2-adic valuations are v2(d_4)=3, v2(d_3)=2, v2(d_2)=2, v2(d_1)=0, v2(d_0)=3.
    # The Newton polygon points (i, v2(d_i)) are (4,3), (3,2), (2,2), (1,0), (0,3).
    # The lower convex hull has two segments:
    # 1. From (1,0) to (0,3): slope = -3. Horizontal length = 1. This means 1 root with valuation -(-3) = 3.
    # 2. From (1,0) to (4,3): slope = 1. Horizontal length = 3. This means 3 roots with valuation -1.
    print("The valuations of the roots of f(x) are:")
    print(" - One root at x=0 (positive valuation).")
    print(" - From the Newton Polygon, we find one root with 2-adic valuation 3 (positive).")
    print(" - From the Newton Polygon, we find three roots with 2-adic valuation -1 (negative).")
    print("-" * 20)
    
    # Step 5: Group branch points into clusters
    # Branch points reducing to 0: roots with positive valuation.
    # Branch points reducing to infinity: roots with negative valuation, plus the point at infinity.
    b0_size = 2  # (root 0, root with v2=3)
    b_inf_size = 3 + 1  # (three roots with v2=-1, plus infinity)
    k = 2  # Number of clusters
    print("We group the branch points into clusters based on their 2-adic reduction:")
    print(f"Cluster B_0 (reducing to 0): {b0_size} branch points.")
    print(f"Cluster B_inf (reducing to infinity): {b_inf_size} branch points.")
    print(f"This gives k = {k} components in the stable reduction.")
    print("-" * 20)
    
    # Step 6: Calculate the genus of each component
    g0 = math.floor((b0_size - 1) / 2)
    g_inf = math.floor((b_inf_size - 1) / 2)
    sum_gi = g0 + g_inf
    print("The geometric genus of each component C_i is calculated as g_i = floor((|B_i|-1)/2):")
    print(f"Genus of C_0: g_0 = floor(({b0_size}-1)/2) = {g0}")
    print(f"Genus of C_inf: g_inf = floor(({b_inf_size}-1)/2) = {g_inf}")
    print("-" * 20)
    
    # Step 7 & 8: Apply the genus formula and solve for the number of double points (delta)
    # Formula: g = sum(g_i) + delta - k + 1
    # delta = g - sum(g_i) + k - 1
    num_double_points = g - sum_gi + k - 1
    
    print("Using the formula for the genus of a reducible curve:")
    print("g = (sum of component genera) + num_double_points - num_components + 1")
    print("We can now solve for the number of double points:")
    print(f"{g} = ({g0} + {g_inf}) + num_double_points - {k} + 1")
    print(f"{g} = {sum_gi} + num_double_points - {k} + 1")

    # Simplify the equation
    rhs_constant = sum_gi - k + 1
    print(f"{g} = num_double_points + {rhs_constant}")
    print(f"num_double_points = {g} - {rhs_constant} = {num_double_points}")
    print("-" * 20)
    print(f"The number of double points is {num_double_points}.")
    
    return num_double_points

if __name__ == '__main__':
    solve_curve_reduction()
