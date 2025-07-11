import math

def solve_trajectory():
    """
    Calculates the particle's position at t = 2*sqrt(3) based on the provided
    initial conditions and the simplified classical-like trajectory.
    """
    
    # Final time given in the problem
    t = 2 * math.sqrt(3)
    
    # The initial position x(0) is given by a complex formula.
    # Let's calculate its components.
    x0_constant_part = 3.0
    
    # The term we denote 'u' in the explanation
    # u = (6*(3 - sqrt(3)))^(1/3) + (6*(3 + sqrt(3)))^(1/3)
    term1_base = 6 * (3 - math.sqrt(3))
    term2_base = 6 * (3 + math.sqrt(3))
    x0_variable_part = term1_base**(1/3) + term2_base**(1/3)

    # Full initial position x(0)
    x0 = x0_constant_part + x0_variable_part

    # As per our analysis, the trajectory follows x(t) = x(0) - t^2 / 4.
    # Let's calculate the displacement term delta_x = t^2 / 4.
    t_squared = t**2
    delta_x = t_squared / 4

    # The final position x(t)
    x_final = x0 - delta_x
    
    # Output the steps as required.
    print(f"The initial position is given by x(0) = 3 + (6*(3 - sqrt(3)))^(1/3) + (6*(3 + sqrt(3)))^(1/3)")
    print(f"Let u = (6*(3 - sqrt(3)))^(1/3) + (6*(3 + sqrt(3)))^(1/3) = {x0_variable_part:.10f}")
    print(f"So, x(0) = 3 + u = {x0:.10f}")
    print("\nAssuming the particle follows a classical-like trajectory with zero initial velocity, x(t) = x(0) - t^2/4.")
    print(f"The final time is t = 2*sqrt(3)")
    print(f"The displacement term is -t^2/4 = -(2*sqrt(3))^2 / 4 = -{t_squared:.2f}/4 = -{delta_x:.2f}")

    print("\nThe final position is x_final = x(0) - 3")
    print(f"The equation for the final position is:")
    # We explicitly show each number in the final equation.
    print(f"x(2*sqrt(3)) = {x0:.10f} - {delta_x:.2f}")

    # The result simplifies to (3 + u) - 3 = u.
    print(f"\nThe calculated final position is: {x_final:.10f}")

solve_trajectory()