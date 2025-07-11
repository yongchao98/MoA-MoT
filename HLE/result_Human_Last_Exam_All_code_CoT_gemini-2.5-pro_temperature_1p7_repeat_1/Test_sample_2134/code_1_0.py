import math

def solve_trajectory():
    """
    Calculates the particle's position x(t) at t = 2*sqrt(3).
    
    The trajectory is assumed to follow x(t) = x(0) + t^2/4, based on the classical
    motion under the external potential.
    """
    
    # Time at which to evaluate the position
    t = 2 * math.sqrt(3)
    
    # Calculate the components of the initial position x(0)
    # x(0) = 3 + (6*(3-sqrt(3)))^(1/3) + (6*(3+sqrt(3)))^(1/3)
    # This can also be written as:
    # x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
    
    sqrt3 = math.sqrt(3)
    
    # Base for the first cubic root term
    term1_base = 18 - 6 * sqrt3
    # Base for the second cubic root term
    term2_base = 18 + 6 * sqrt3
    
    # Calculate the cubic root terms
    cbrt_term1 = term1_base**(1/3)
    cbrt_term2 = term2_base**(1/3)
    
    # Calculate the initial position x(0)
    x0_val = 3 + cbrt_term1 + cbrt_term2
    
    # Calculate the displacement term t^2 / 4
    displacement = t**2 / 4
    
    # Calculate the final position x(t)
    xt_val = x0_val + displacement
    
    # Print the breakdown of the final equation as requested
    print(f"The trajectory is determined by the equation: x(t) = x(0) + t^2 / 4")
    print(f"The initial position is x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)")
    print(f"Numerically, x(0) = {x0_val}")
    print(f"At t = 2*sqrt(3), the displacement term is t^2 / 4 = {displacement}")
    print(f"The final position is x(2*sqrt(3)) = x(0) + 3 = {x0_val} + {displacement}")
    print(f"x(2*sqrt(3)) = {xt_val}")

solve_trajectory()