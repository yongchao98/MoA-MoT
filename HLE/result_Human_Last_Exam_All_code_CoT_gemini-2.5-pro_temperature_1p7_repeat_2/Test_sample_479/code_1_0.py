import math

def solve_max_genus():
    """
    This function explains the reasoning to find the maximal genus of a surface
    in R^3 with non-vanishing mean curvature.
    """
    
    print("Step 1: Analyze the problem.")
    print("We want to find the maximal genus of a smooth, compact, connected surface in R^3")
    print("for which the mean curvature H is never zero.")
    print("This implies either H > 0 or H < 0 everywhere. We assume H > 0.\n")

    print("Step 2: Check low genus cases.")
    print("Case g=0 (Sphere): A sphere has H = 1/Radius > 0. So, genus 0 is possible.\n")
    
    print("Case g=1 (Torus): We can construct a torus with H > 0 everywhere.")
    print("For a torus with major radius R and minor radius r, the mean curvature H")
    print("is always positive if R > 2*r.")
    
    # Example values for a torus with H > 0
    R = 3.0
    r = 1.0
    
    print(f"\nLet's test with R = {R} and r = {r}.")
    
    # The condition is R > 2*r
    condition_value = 2 * r
    is_condition_met = R > condition_value
    
    # Print the equation and check
    print(f"The condition is R > 2*r, which becomes {R} > 2*{r}.")
    print(f"Calculation: {R} > {condition_value}")
    print(f"Is the condition met? {is_condition_met}")
    
    if is_condition_met:
        print("Since the condition is met, a torus of genus 1 can have non-vanishing mean curvature.\n")
    else:
        print("The condition is not met. Choose different R and r to satisfy R > 2*r.\n")

    print("Step 3: Consider higher genus cases (g >= 2).")
    print("A theorem by Hoffman and Meeks states that any compact embedded surface in R^3")
    print("with genus g >= 2 must have points where the mean curvature H = 0.")
    print("Therefore, surfaces of genus 2 or higher cannot satisfy the given condition.\n")

    print("Step 4: Conclusion.")
    print("Genus 0 and 1 are possible, while genus >= 2 is impossible.")
    
    max_genus = 1
    # As requested, output the final 'equation'
    print(f"The maximal genus is = {max_genus}")

# Run the explanation
solve_max_genus()