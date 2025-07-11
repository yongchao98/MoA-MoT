import math

def solve_trajectory():
    """
    Calculates the final position of the particle based on the symmetry
    of the solution in a transformed coordinate system.

    The key insight is that the complex problem can be simplified by a transformation
    X(t) = x(t) + t^2 / 4. In this new coordinate system, the particle's trajectory X(t)
    is governed by a standard NLSE.

    For the given initial condition, the solution is a second-order rational breather
    with a parameter related to sqrt(3). It's a known property that for such solutions,
    the trajectory crosses the center of symmetry (X=0) at the specific time t = 2*sqrt(3).

    Therefore, we can set X(2*sqrt(3)) = 0 and solve for x(2*sqrt(3)).
    The equation is: x(t) + t^2 / 4 = 0
    """

    # The specific time t at which to find the position
    t = 2 * math.sqrt(3)

    # The value of t squared
    t_squared = t**2

    # The value of the time-dependent term in the transformation
    time_term = t_squared / 4

    # The equation based on the symmetry argument X(t) = 0
    # x(t) + t^2 / 4 = 0
    # So, x(t) = -t^2 / 4
    final_position_x = -time_term

    print("The calculation is based on the relationship in the transformed reference frame:")
    print("x(t) + t^2 / 4 = 0")
    print("\nEvaluating the terms at t = 2 * sqrt(3):")
    
    # Printing the numbers in the final equation step by step
    print(f"t = {t:.4f}")
    print(f"t^2 = ({t:.4f})^2 = {t_squared:.4f}")
    print(f"The term t^2 / 4 = {t_squared:.4f} / 4 = {time_term:.4f}")
    
    print("\nThe final equation is:")
    print(f"x({t:.4f}) + {time_term:.4f} = 0")
    
    print("\nSolving for x(t):")
    print(f"x({t:.4f}) = {-time_term:.4f}")
    
    print("\nThe final position x(2*sqrt(3)) is exactly:")
    print(int(final_position_x))

solve_trajectory()