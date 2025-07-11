import math

def explain_maximal_genus():
    """
    Explains the reasoning for determining the maximal genus of a specific type of surface.
    """
    print("Step 1: Determine the sign of the Mean Curvature (H)")
    print("A standard result from differential geometry (the height function argument) shows that for any compact surface bounding a region in R^3,")
    print("there must be at least one point 'p_max' (a maximum along some direction) where H(p_max) <= 0.")
    print("The problem states that the mean curvature vector never vanishes, so H is never 0.")
    print("Since the surface is connected, H must have a constant sign everywhere. Combined with H(p)<=0 at some point, this means H must be strictly negative (H < 0) everywhere.")
    print("\nStep 2: Check if Genus 1 is possible")
    print("Let's consider a torus (genus 1) as the boundary of a solid torus (a compact region).")
    print("Let the major radius be R and the minor radius be r.")
    print("The mean curvature H (with respect to the outward normal) is given by the formula:")
    print("H(theta) = - (R + 2*r*cos(theta)) / (2*r*(R + r*cos(theta)))")
    print("where theta is the poloidal angle.")

    print("\nTo have H < 0 everywhere, the numerator must be positive everywhere, since the denominator is always positive for a torus (R > r).")
    print("We need: R + 2*r*cos(theta) > 0 for all theta.")
    print("The minimum value of this expression occurs at cos(theta) = -1.")
    print("So we require: R + 2*r*(-1) > 0  => R - 2*r > 0 => R > 2*r.")

    R = 3
    r = 1
    print(f"\nFor example, let's choose R = {R} and r = {r}. This satisfies the condition R > 2*r ({R} > 2*{r}).")
    
    # Calculate H at different points
    theta_out = 0  # outside
    H_out = - (R + 2*r*math.cos(theta_out)) / (2*r*(R + r*math.cos(theta_out)))
    print(f"At the outermost point (theta=0), H = -({R} + 2*{r}*1) / (2*{r}*({R} + {r}*1)) = {H_out:.4f}, which is < 0.")

    theta_in = math.pi  # inside
    H_in = - (R + 2*r*math.cos(theta_in)) / (2*r*(R + r*math.cos(theta_in)))
    print(f"At the innermost point (theta=pi), H = -({R} + 2*{r}*(-1)) / (2*{r}*({R} + {r}*(-1))) = {H_in:.4f}, which is < 0.")
    
    print("\nSince we can construct a valid genus 1 surface, the maximal genus is at least 1.")

    print("\nStep 3: Higher Genera (g >= 2)")
    print("A deeper result in geometry shows that it's impossible for an embedded surface of genus g >= 2 that bounds a compact region in R^3 to have its mean curvature be of a single sign (everywhere positive or everywhere negative).")
    print("Such surfaces are topologically more complex and are forced to have H=0 at some points.")
    
    print("\nConclusion:")
    print("The maximal genus is 1.")

if __name__ == '__main__':
    explain_maximal_genus()
